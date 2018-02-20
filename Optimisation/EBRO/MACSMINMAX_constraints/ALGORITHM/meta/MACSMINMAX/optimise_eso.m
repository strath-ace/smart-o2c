% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [d,fval,exitflag,output] = optimise_eso(problem_minmax, algo_outer, algo_inner, par_minmax)

global nfevalglobal
% Rename general inputs
n_d = problem_minmax.dim_d;
n_u = problem_minmax.dim_u;
n_obj = problem_minmax.n_obj;
sign_inner = problem_minmax.sign_inner; % 1 for minmax, -1 for minmin
nfevalmax = par_minmax.maxnfeval;

% Rename MinMaReK inputs
indicator_d_max = par_minmax.indicator_d_max; % Threshold on the indicator (EI, PI, etc.) for outer loop
indicator_u_max = par_minmax.indicator_u_max; % Threshold on the indicator (EI, PI, etc.) for inner loop
iter_d_max = par_minmax.iter_d_max; % Max. iterations for outer loop in an iteration of MinMaReK
iter_u_max = par_minmax.iter_u_max; % Max. iterations for inner loop in an iteration of MinMaReK
tol_conv = par_minmax.tol_conv; % Tolerance to convergence, only implemented for S.O.

lsflag_validation = par_minmax.local_search_flags.validation; % run either local search or func eval in validation
keep_d_record = par_minmax.keep_d_record;

% Retrieve initial surrogate
surrogate = par_minmax.surrogate;              % CHANGED init_algo_emo

% Hardcoded MinMaReK inputs
n_du_0 = 10*(n_d+n_u); %HC
n_u_record_0 = repmat([1],[1,n_obj]); %HC

% Build a mock macsminmax metaproblem max_u that we will use for validation and function evaluation
problem_max_u = build_metaproblem_macsminmax_inner(problem_minmax);

% Build true inner and outer problems of MinMaReK
% This will also initialise the surrogates
problem_outer = build_metaproblem_ideaminmax_s_outer(problem_minmax);
problem_inner = build_metaproblem_ideaminmax_s_inner(problem_minmax);
problem_outer.par_objfun.indicator_d = par_minmax.indicator_d;

% The fun begins
nfeval = 0;

% Initialise DOE
du_0 = lhsu(zeros(1,n_d+n_u),ones(1,n_d+n_u),n_du_0);
x_doe = [];
f_doe = [];
for i=1:n_du_0
    x_doe(i,:) = du_0(i,:);
    problem_max_u.par_objfun.d = du_0(i,1:n_d);
    % NOTE: A lot of the fevals will look like this because so far the code takes one 
    % function handle per objective. Some adaptation is needed so it can take vectorial objfuns.
    for obj2 = 1:n_obj 
        problem_max_u.par_objfun.objective = obj2;
        f_doe(i,obj2) = -sign_inner*problem_max_u.objfun(du_0(i,n_d+1:n_d+n_u),problem_max_u.par_objfun);
    end
    nfeval = nfeval + 1; % Conceptually this is only one feval.
end

% train surrogate inner

% problem_inner.par_objfun.surrogate.model = problem_inner.par_objfun.surrogate.training(x_doe,-sign_inner*f_doe,problem_inner.par_objfun.surrogate);
for obj = 1:n_obj
    surrogate{obj} = surrogate{obj}.add(x_doe, f_doe(:,obj), surrogate{obj});
    surrogate{obj} = surrogate{obj}.update(surrogate{obj});
end
problem_inner.par_objfun.surrogate = surrogate;

% Declare the optimisation archives Au, Ad, Afminmax
u_record = cell(1,n_obj);
d_record = [];
f_record = [];

% % Initialise archive Au randomly
% for obj = 1:n_obj
%     u_record{obj} = lhsu(zeros(1,n_u),ones(1,n_u),n_u_record_0(obj));
% end

% Initialise archive Au as they do it in ideaminmax_sur3:
for obj = 1:n_obj
    % look for the d with max f in the DOE.
    % NOTE: what if we try with the min??
    [~, idx] = max(f_doe(:,obj));
    dmaxdoe = x_doe(idx,1:n_d);
    
    % run an EI optimisation in U for those ds
    problem_inner.par_objfun.d = dmaxdoe;
    problem_inner.par_objfun.objective = obj;
    problem_inner.par_objfun.ymin = 1e32*ones(1,n_obj);

    [umax_dmaxdoe , ~ , ~ , ~ ] = algo_inner.optimise(problem_inner,algo_inner.par);

    % add the output u to Au
    u_record{obj}=[u_record{obj};umax_dmaxdoe];

    % evaluate the output x and add to the DOE
    f_aux = [];
    problem_max_u.par_objfun.d = dmaxdoe;
    for obj2 = 1:n_obj
        problem_max_u.par_objfun.objective = obj2;
        f_aux(1,obj2) = -sign_inner*problem_max_u.objfun(umax_dmaxdoe,problem_max_u.par_objfun);
        surrogate{obj} = surrogate{obj}.add([dmaxdoe umax_dmaxdoe], f_aux(1,obj2), surrogate{obj});
    end
    nfeval = nfeval + 1;
    
    x_doe = [x_doe;[dmaxdoe umax_dmaxdoe]];
    f_doe = [f_doe;f_aux];

end



%% Main loop
iter=0;
precision = inf;
dmin_outer = nan(1,n_d);
fmin_outer = 1e32*ones(1,n_obj);
umax_outer=cell(1,n_obj);
umax_inner_history = cell(1,n_obj);
for obj = 1:n_obj
    umax_outer{obj}=nan(1,n_u);
    umax_inner_history{obj}= nan(1,n_u);
end
fmax_inner_history = 1e32*ones(1,n_obj);
stop = false;
while ~stop

    iter = iter +1;
    u_record_aux = cell(1,n_obj);

    %% OUTER LOOP: MINIMISATION OVER D
    indicator_d = inf;
    iter_d = 0;
    problem_outer.par_objfun.u_record = u_record;
    while abs(indicator_d) > indicator_d_max && iter_d < iter_d_max
        
        iter_d = iter_d+1;

        % update surrogate of outer problem
        [~, idx] = unique(round(1e8*x_doe),'rows');
        x_doe = x_doe(idx,:);
        f_doe = f_doe(idx,:);
      
        % problem_outer.par_objfun.surrogate.model = problem_outer.par_objfun.surrogate.training(x_doe,f_doe,problem_outer.par_objfun.surrogate);
        for obj = 1:n_obj
            surrogate{obj} = surrogate{obj}.update(surrogate{obj});
        end
        problem_outer.par_objfun.surrogate = surrogate;

        problem_outer.par_objfun.ymin = fmin_outer; % we'll come back to this, but it's like this in ideaminmax_sur3

        % maximise indicator_d
        [d_outer_aux, indicator_d, ~, ~] = algo_outer.optimise(problem_outer,algo_outer.par); %no feval here
        indicator_d = -indicator_d;

        % update x_doe, f_doe (1)
        % also build f_outer_aux that is like the result of validation on d_outer_aux and u_record
        % (running u_validation would repeat fevals). This indeed looks like an inline u_validation.
        x_doe_aux = [];
        f_doe_aux = [];
        f_outer_aux = -sign_inner*inf(size(d_outer_aux,1),n_obj);
        u_outer_aux = cell(1,n_obj);
        for obj = 1:n_obj
            u_outer_aux{obj} = nan(size(d_outer_aux,1),problem_minmax.dim_u);
        end
        % We will add dmin x Au to the DOE
        for idx_d = 1:size(d_outer_aux,1) % this should be 1:1 but let's leave it like this to facilitate possible changes
            problem_max_u.par_objfun.d = d_outer_aux(idx_d,:);
            for obj = 1:n_obj
                for idx_u = 1:size(u_record{obj})
                    f_doe_aux_aux =[];
                    for obj2 = 1:n_obj
                        problem_max_u.par_objfun.objective = obj2;
                        f_doe_aux_aux(1,obj2) = -sign_inner*problem_max_u.objfun(u_record{obj}(idx_u,:),problem_max_u.par_objfun);
                    end
                    nfeval = nfeval + 1;

                    x_doe_aux = [x_doe_aux; d_outer_aux(idx_d,:), u_record{obj}(idx_u,:)];
                    f_doe_aux = [f_doe_aux; f_doe_aux_aux];

                    % And in the meanwhile we keep track of the max. in Au
                    if (sign_inner*f_outer_aux(idx_d,obj)<sign_inner*f_doe_aux_aux(1,obj))
                        f_outer_aux(idx_d,obj) = f_doe_aux_aux(1,obj);
                        u_outer_aux{obj}(idx_d,:) = u_record{obj}(idx_u,:);
                    end
                end
            end
        end
        x_doe = [x_doe; x_doe_aux];
        f_doe = [f_doe; f_doe_aux];
        for obj = 1:n_obj
            surrogate{obj} = surrogate{obj}.add(x_doe_aux, f_doe_aux(:,obj), surrogate{obj});
        end

        % update fmin_outer for EI computation, the idea is that in the MO case a minmax front gets gradually built
        idx=[];
        if n_obj == 1
            [fmin_outer,idx] = min([fmin_outer;f_outer_aux]);
            dmin_outer_stack = [dmin_outer; d_outer_aux];
            dmin_outer = dmin_outer_stack(idx,:);
        else
            f_outer_stack = [fmin_outer; f_outer_aux];
            d_outer_stack = [dmin_outer; d_outer_aux];
            idx = dominance(f_outer_stack,0) == 0;
            fmin_outer = f_outer_stack(idx,:);
            dmin_outer = d_outer_stack(idx,:);
        end

        for obj=1:n_obj
            umax_outer_stack_obj = [umax_outer{obj};u_outer_aux{obj}];
            umax_outer{obj} = umax_outer_stack_obj(idx,:);
        end


        % OLD VERSION OF WHAT IS ABOVE SINCE (1)
        % % validate dmin with the archive
        % [f_outer_aux, u_outer_aux, nfeval_aux, all_f] = u_validation(problem_max_u, d_outer_aux, u_record, false, 1:n_obj);
        % nfeval = nfeval + nfeval_aux;

        % % update x_doe, f_doe and surrogate
        % % NOTE suppose there is one doe per objective and we'll see if this works or not...
        % % NOTE II in fact this is stupid because evaluating one objective should be more or less as hard as evaluating them all... 
        % % so I should recode with just one doe asap ("recoding" validation here) but for the moment leave it like this for SO version
        % for obj = 1:n_obj
        %     size_u_record_obj = size(u_record{obj},1);
        %     for i = 1:size(d_outer_aux,1)
        %         x_doe{obj} = [x_doe{obj}; repmat(d_outer_aux(i,:),[size_u_record_obj,1]), u_record{obj}];
        %         f_doe{obj} = [f_doe{obj}; all_f{obj}(i,:)'];
        %     end
        %     % [x_doe{obj},idx] = uniquetol(x_doe{obj},1e-8,'ByRows',true); % R2015+
        %     [~, idx] = unique(round(1e8*x_doe{obj}),'rows');
        %     x_doe{obj} = x_doe{obj}(idx,:);
        %     f_doe{obj} = f_doe{obj}(idx,:);
        %     problem_outer.par_objfun.surrogate.model{obj} = problem_outer.par_objfun.surrogate.training(x_doe{obj},f_doe{obj},problem_outer.par_objfun.surrogate);
        % end

        % if n_obj == 1
        %     [fmin_outer,idx] = min([fmin_outer;f_outer_aux]);
        %     dmin_outer = [dmin_outer; d_outer_aux];
        %     dmin_outer = dmin_outer(idx,:);
        % else
        %     error('MO: some things not implemented! (and nothing tested...)')
        %     % problem_outer.par_objfun.ymin = f_outer(dominance(f_outer,0) == 0,:); % pareto front
        %     % probably the simpler is to make a stack f_outer with all f_outer_aux and choose the non-dominated
        % end
        
    end

    % fmin_outer

    % % Archive shrinking for MO??
    
    %% INNER LOOP: MAXIMISATION OVER U
    n_dmin = size(dmin_outer,1);
    fmax_inner = fmin_outer;
    umax_inner = umax_outer;
    f_record_aux=fmin_outer;

    for i = 1:n_dmin
        problem_inner.par_objfun.d = dmin_outer(i,:);
        problem_max_u.par_objfun.d = dmin_outer(i,:);
        list_obj = randperm(n_obj);
        for obj = list_obj
            problem_max_u.par_objfun.objective = obj;
            problem_inner.par_objfun.objective = obj;
            % fmax_inner_obj = -sign_inner*fmin_outer(i,obj); % scalar! % as it should be according to MinMaReK
            % umax_inner_obj = umax_outer{obj}(i,:);
            fmax_inner_obj = fmax_inner_history(1,obj); % scalar!
            umax_inner_obj = umax_inner_history{obj};
            indicator_u = inf;
            iter_u = 0;
            while abs(indicator_u) > indicator_u_max && iter_u < iter_u_max
                iter_u = iter_u + 1;
                % train surrogate of inner problem
                % NOTE I  : the 2nd objective has more info than the 1st, etc... potential issue?
                % NOTE II : the inner surrogate gets trained in negative for minmax. maybe it should be trained in positive
                %           and just the EI computed for maximisation, so the surrogate can be passed from inner to outer
                [~, idx] = unique(round(1e8*x_doe),'rows');
                x_doe = x_doe(idx,:);
                f_doe = f_doe(idx,:);

%                 for obj2 = 1:n_obj
                surrogate{obj} = surrogate{obj}.update(surrogate{obj});
%                 end

                %problem_inner.par_objfun.surrogate.model = problem_inner.par_objfun.surrogate.training(x_doe,-sign_inner*f_doe,problem_inner.par_objfun.surrogate);
                problem_inner.par_objfun.surrogate = surrogate;
                problem_inner.par_objfun.ymin = fmax_inner_obj*ones(1,n_obj); % we'll come back to this, but it's like this in ideaminmax_sur3

                % maximise indicator_u
                [u_inner_aux, indicator_u, ~, ~] = algo_inner.optimise(problem_inner,algo_inner.par); %no nfeval here
                indicator_u = -indicator_u;

                % Add the point to the DOE
                f_doe_aux = [];
                for obj2 = 1:n_obj
                    f_doe_aux(1,obj2) = problem_max_u.objfun(u_inner_aux,problem_max_u.par_objfun); %note it comes negative for minmax
                    surrogate{obj2} = surrogate{obj2}.add([dmin_outer(i,:) , u_inner_aux], -sign_inner*f_doe_aux(1,obj2), surrogate{obj2});
                end
                nfeval = nfeval + 1;

                f_inner_aux_obj = f_doe_aux(1,obj); % scalar!
                
                % update the DOE with the new point
                x_doe = [x_doe; dmin_outer(i,:) , u_inner_aux];
                f_doe = [f_doe; -sign_inner*f_doe_aux];

                % update ymin for EI computation
                [fmax_inner_obj,idx] = min([fmax_inner_obj; f_inner_aux_obj]); % scalar!
                umax_inner_obj_stack = [umax_inner_obj; u_inner_aux];
                umax_inner_obj = umax_inner_obj_stack(idx,:);
            end
            
            fmax_inner(i,obj) = -sign_inner*fmax_inner_obj;
            
            % update Au with the champion between outer and inner loop
            if fmax_inner_obj < -sign_inner*fmin_outer(i,obj)
                u_record_aux{obj} = [u_record_aux{obj}; umax_inner_obj];
                umax_inner{obj}(i,:)=umax_inner_obj;
                f_record_aux(i,obj) = -sign_inner*fmax_inner_obj;
            else
                if iter==1 %this is already in the archive, only replace it for iter==1 because we are gonna restart the archive
                    u_record_aux{obj} = [u_record_aux{obj}; u_outer_aux{obj}];
                end
            end
        end
    end
    d_record = [d_record;dmin_outer];
    f_record = [f_record;f_record_aux];
    
    % f_record_aux

    % update u_record
    for obj = 1:n_obj
        if iter==1
            u_record{obj} = []; % NOTE: not sure if this is good for the MO but was like this in ideaminmax_sur3
        end
        u_record{obj} = [u_record{obj}; u_record_aux{obj}];
        % u_record{obj} = uniquetol(u_record{obj},1e-8,'ByRows',true); % R2015b+
        [~,idx] = unique(round(1e8*u_record{obj}),'rows');
        u_record{obj} = u_record{obj}(idx,:);
    end
    
    % Stop criterions
    size_u_record = sum(cellfun('size',u_record,1));
    nfeval_val = size_u_record*size(dmin_outer,1); % NOTE: Should add a rule of thumb.
    stop = nfeval + nfeval_val >= nfevalmax;
    if n_obj == 1
        % NOTE: Maybe the idea of a precision in the maximisation as stop criterion can be extended to MO?
        precision = abs(-sign_inner*fmax_inner_obj - fmin_outer);
        stop = stop || precision < tol_conv;
    end 

    % update the f-values that we use to compute EI
    % NOTE I  ; Needs to be done here not to interfer with stopping criterion for SO
    % NOTE II : Reinitialised the opposite to add randomness to the iteration?? Should test
    %           the recommended initialisation of MinMaReK maybe improves convergence, 
    %           especially in MO but for the moment I follow ideaminmax_sur3;
    [fmax_inner_history,idx] = min(sign_inner*f_doe,[],1);
    fmax_inner_history = -sign_inner*fmax_inner_history;
    for obj = 1:n_obj
        umax_inner_history{obj} = x_doe(idx(obj),n_d+1:n_d+n_u);
    end
    idx=[];
    if n_obj == 1
        [fmin_outer,idx] = max(f_doe,[],1);
        if(~stop)
            dmin_outer = x_doe(idx,1:n_d);
        end
        for obj = 1:n_obj
            umax_outer{obj} = x_doe(idx,n_d+1:n_d+n_u);
        end
    else
        % % no modification on dmin outer
        % fmin_outer = f_record_aux;
        % umax_outer = umax_inner;

        idx = dominance(-f_doe,0)==0;
        fmin_outer = f_doe(idx,:);
        if(~stop); dmin_outer = x_doe(idx,1:n_d); end;
        for obj = 1:n_obj
            umax_outer{obj} = x_doe(idx,n_d+1:n_d+n_u);
        end
    end
    
    
    % iter
    % nfeval
    % size_x_doe = size(x_doe,1)

end

if ~keep_d_record
    d_record = dmin_outer;
end

%% Archive cross-check
% NOTE: even more clever would be to keep track of which u's have the d's been validated against already
stop = false;
checked = false(1,size(d_record,1));
sel = false(1,size(d_record,1));
u_val_record = cell(1,n_obj);
for obj = 1:n_obj
    u_val_record{obj}=nan(size(d_record,1),problem_minmax.dim_u);
end
nfeval_val = 0;
while ~stop
    sel = dominance(f_record,0) == 0;       % find non-dominated
    if(n_obj>1) 
        sel = sel';
    end
    
    tocheck = sel & ~checked;               % select only those that have not been checked
    stop = all(~tocheck);                   % stop if you have nothing to check
    
    d_tocheck = d_record(tocheck,:);        % check those that have been selected

    % validate
    [f_val_aux, u_val_record_aux, nfeval_aux] = u_validation(problem_max_u, d_tocheck, u_record, lsflag_validation, 1:n_obj);
    nfeval_val = nfeval_val + nfeval_aux;

    % update f_record and u_val_record
    f_record(tocheck,:) = f_val_aux;
    for obj = 1:n_obj
        u_val_record{obj}(tocheck,:) = u_val_record_aux{obj};
    end
    checked = checked | tocheck;           % update who has been checked
end

% at the end of the loop above, sel keeps the non-dominated and validated solutions
fval = f_record(sel,:);
nfeval = nfeval + nfeval_val;

%%
frontsize = size(fval,1);

d = d_record(sel,:).*repmat(problem_minmax.ub_d'-problem_minmax.lb_d',[frontsize,1]) + repmat(problem_minmax.lb_d',[frontsize,1]);
u = cell(1,n_obj);
for obj = 1:n_obj
    u_val_record{obj} = u_val_record{obj}(sel,:);
    map_info = problem_max_u.par_objfun.map_u_info{obj};
    for i = 1:frontsize
        u{obj}(i,:) = map_affine(u_val_record{obj}(i,:),map_info); %this can be easily vectorized
    end
end

output.u = u;
output.nfeval = nfeval;
exitflag = 0;
if n_obj == 1
    exitflag = nfeval <= nfevalmax && precision < tol_conv; 
end

end

