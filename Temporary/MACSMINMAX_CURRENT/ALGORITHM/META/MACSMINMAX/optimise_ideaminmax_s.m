function [d,fval,exitflag,output] = optimise_ideaminmax_s(problem_minmax, algo_outer, algo_inner, par_minmax)

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

% Hardcoded MinMaReK inputs
n_du_0 = 10*(n_d+n_u); %HC
n_u_record_0 = repmat([1],[1,n_obj]); %HC

% Build a mock macsminmax metaproblem max_u that we will use for validation and function evaluation
problem_max_u = build_metaproblem_macsminmax_inner(problem_minmax);

% Build true inner and outer problems of MinMaReK
% This will also initialise the surrogates
problem_outer = build_metaproblem_ideaminmax_s_outer(problem_minmax);
problem_inner = build_metaproblem_ideaminmax_s_inner(problem_minmax);

% The fun begins
nfeval = 0;

% Initialise DOE and train surrogate
du_0 = lhsu(zeros(1,n_d+n_u),ones(1,n_d+n_u),n_du_0);
x_doe = [];
f_doe = [];
for i=1:n_du_0
    x_doe(i,:) = du_0(i,:);
    problem_max_u.par_objfun.d = du_0(i,1:n_d);
    for obj = 1:n_obj % This is like this because so far the code takes one function handle per objective. Some adaptations are needed so it can take vectorial objfuns.
        problem_max_u.par_objfun.objective = obj;
        f_doe(i,obj) = -sign_inner*problem_max_u.objfun(x_doe(i,n_d+1:n_d+n_u),problem_max_u.par_objfun);
    end
    nfeval = nfeval + 1;
end

% train surrogate inner
problem_inner.par_objfun.surrogate.model = problem_inner.par_objfun.surrogate.training(x_doe,-sign_inner*f_doe,problem_inner.par_objfun.surrogate);


u_record = cell(1,n_obj);
d_record = [];
f_record = [];

% % Initialise archive Au randomly and Ad empty
% % NOTE: instead of this in ideaminmax here they select the minimum point in x_doe and take their d, 
% % optimise in u, add us to the archive Au, ds and us to the x_doe and retrain.
% for obj = 1:n_obj
%     u_record{obj} = lhsu(zeros(1,n_u),ones(1,n_u),n_u_record_0(obj));
% end

% as they do it in ideaminmax_sur3:
problem_inner.par_objfun.ymin = 1e32;
for obj = 1:n_obj
    
    [~, idmaxdoe] = max(f_doe(:,obj));
    dmaxdoe = x_doe(idmaxdoe,1:n_d);
    
    problem_inner.par_objfun.d = dmaxdoe;
    problem_inner.par_objfun.objective = obj;
    
    [umaxmaxdoe , ~ , ~ , ~ ] = algo_inner.optimise(problem_inner,algo_inner.par);

    f_aux = [];
    problem_max_u.par_objfun.d = dmaxdoe;
    for obj2 = 1:n_obj
        problem_max_u.par_objfun.objective = obj2;
        f_aux(1,obj2) = -sign_inner*problem_max_u.objfun(umaxmaxdoe,problem_max_u.par_objfun);
    end
    nfeval = nfeval + 1;
    
    u_record{obj}=[u_record{obj};umaxmaxdoe];
    
    x_doe = [x_doe;[dmaxdoe umaxmaxdoe]];
    f_doe = [f_doe;f_aux];
end



% Main loop
iter=0;
precision = inf;
problem_outer.par_objfun.ymin = 1e32;
% problem_inner.par_objfun.ymin = 1e32;
dmin_outer = nan(1,n_d);
fmin_outer = 1e32*ones(1,n_obj);
umax_inner = nan(1,n_d);
fmax_inner = sign_inner*1e32*ones(1,n_obj);
stop = false;
while ~stop

    iter = iter +1;
    u_record_aux = cell(1,n_obj);

    % %OUTER LOOP: MINIMISATION OVER D
    indicator_d = inf;
    iter_d = 0;
    problem_outer.par_objfun.u_record = u_record;
    while abs(indicator_d) > indicator_d_max && iter_d < iter_d_max

        iter_d = iter_d+1;

        % update surrogate of outer problem
        [~, idx] = unique(round(1e8*x_doe),'rows');
        x_doe = x_doe(idx,:);
        f_doe = f_doe(idx,:);
        problem_outer.par_objfun.surrogate.model = problem_outer.par_objfun.surrogate.training(x_doe,f_doe,problem_outer.par_objfun.surrogate);
        problem_outer.par_objfun.ymin = fmin_outer; % we'll come back to this, but it's like this in ideaminmax

        % maximise indicator_d
        [d_outer_aux, indicator_d, ~, ~] = algo_outer.optimise(problem_outer,algo_outer.par); %no feval here
        indicator_d = -indicator_d;

        % update x_doe, f_doe
        % also build f_outer_aux that is like the result of validation on d_outer_aux and u_record
        % (running validation would repeat fevals). This is indeed like an inline u_validation.
        x_doe_aux = [];
        f_doe_aux = [];
        f_outer_aux = -sign_inner*inf(size(d_outer_aux,1),n_obj);
        u_outer_aux = cell(1,n_obj);
        for obj = 1:n_obj
            u_outer_aux{obj} = nan(size(d_outer_aux,1),problem_minmax.dim_u);
        end

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

                    if (sign_inner*f_outer_aux(idx_d,obj)<sign_inner*f_doe_aux_aux(1,obj))
                        f_outer_aux(idx_d,obj) = f_doe_aux_aux(1,obj);
                        u_outer_aux{obj}(idx_d,:) = u_record{obj}(idx_u,:);
                    end

                    x_doe_aux = [x_doe_aux; d_outer_aux(idx_d,:), u_record{obj}(idx_u,:)];
                    f_doe_aux = [f_doe_aux; f_doe_aux_aux];
                end
            end
        end
        x_doe = [x_doe; x_doe_aux];
        f_doe = [f_doe; f_doe_aux];



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

        % update fmin_outer and friends, and ymin for EI computation
        f_outer_stack = [fmin_outer; f_outer_aux];
        d_outer_stack = [dmin_outer; d_outer_aux];
        idx = dominance(f_outer_stack,0) == 0;
        fmin_outer = f_outer_stack(idx,:);
        dmin_outer = d_outer_stack(idx,:);

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

    n_dmin = size(dmin_outer,1);
    % % Archive shrinking for MO??
    
    % INNER LOOP: MAXIMISATION OVER U
    for i = 1:n_dmin
        problem_inner.par_objfun.d = dmin_outer(i,:);
        problem_max_u.par_objfun.d = dmin_outer(i,:);
        f_record_aux=zeros(1,n_obj);
        for obj = 1:n_obj
            problem_max_u.par_objfun.objective = obj;
            problem_inner.par_objfun.objective = obj;

            indicator_u = inf;
            iter_u = 0;
            while abs(indicator_u) > indicator_u_max && iter_u < iter_u_max
                iter_u = iter_u + 1;

                % train surrogate inner problem
                [~, idx] = unique(round(1e8*x_doe),'rows');
                x_doe = x_doe(idx,:);
                f_doe = f_doe(idx,:);
                problem_inner.par_objfun.surrogate.model = problem_inner.par_objfun.surrogate.training(x_doe,-sign_inner*f_doe,problem_inner.par_objfun.surrogate);
                problem_inner.par_objfun.ymin = fmax_inner(1,obj); % we'll come back to this, but it's like this in ideaminmax

                % maximise indicator_u
                [u_inner_aux, indicator_u, ~, ~] = algo_inner.optimise(problem_inner,algo_inner.par); %no nfeval here
                indicator_u = -indicator_u;

                % evaluate
                f_doe_aux = [];
                for obj2 = 1:n_obj
                    f_doe_aux(1,obj2) = problem_max_u.objfun(u_inner_aux,problem_max_u.par_objfun); %note it comes negative for minmax
                end
                nfeval = nfeval + 1;

                f_inner_aux = f_doe_aux(1,obj);
                
                x_doe = [x_doe; dmin_outer(i,:) , u_inner_aux];
                f_doe = [f_doe; -sign_inner*f_doe_aux];

                % update ymin for EI computation
                if n_obj == 1
                    [fmax_inner,idx] = min([fmax_inner;f_inner_aux]);
                    umax_inner = [umax_inner; u_inner_aux];
                    umax_inner = umax_inner(idx,:);
                else
                    error('MO: some things not implemented! (and nothing tested...)')
                    % problem_outer.par_objfun.ymin = f_outer(dominance(f_outer,0) == 0,:); % pareto front
                end
            end

            if fmax_inner(1,obj) < -sign_inner*fmin_outer(i,obj)
                u_record_aux{obj} = [u_record_aux{obj}; umax_inner];
                f_record_aux(1,obj) = -sign_inner*fmax_inner(1,obj);
            else
                if iter==1 %this is already in the archive, only replace it for iter==1 because we are gonna restart the archive
                    u_record_aux{obj} = [u_record_aux{obj}; u_outer_aux];
                end
                f_record_aux(1,obj) = fmin_outer(i,obj);
            end
        end
        d_record = [d_record;dmin_outer(i,:)];
        f_record = [f_record;f_record_aux];
    end
    
    
    % update u_record
    for obj = 1:n_obj
        if iter==1
            u_record{obj} = []; % I don't like this too much
        end
        u_record{obj} = [u_record{obj}; u_record_aux{obj}];
        % u_record{obj} = uniquetol(u_record{obj},1e-8,'ByRows',true); % R2015b+
        [~,idx] = unique(round(1e8*u_record{obj}),'rows');
        u_record{obj} = u_record{obj}(idx,:);
    end
    
    size_u_record = sum(cellfun('size',u_record,1));
    nfeval_val = size_u_record*size(dmin_outer,1);
    stop = nfeval + nfeval_val >= nfevalmax;
    if n_obj == 1
        precision = abs(-sign_inner*fmax_inner - fmin_outer);
        stop = stop || precision < tol_conv;
    end
    
    fmax_inner = -sign_inner*min(sign_inner*f_doe,[],1);
    fmin_outer = max(f_doe,[],1);
    
end

if ~keep_d_record
    d_record = dmin_outer;
end

%% archive cross check
% [f_val, u_val_record, nfeval_aux] = u_validation(problem_max_u, d_record, u_record, lsflag_validation, 1:n_obj);
% nfeval = nfeval + nfeval_aux;
% % Select non-dominated solutions and define outputs
% dom = dominance(f_val,0) == 0;
% fval = f_val(sel,:);

%% more clever archive cross check
% even more clever would be to keep track of which u's have the d's been validated against already

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
    % check those that have been selected
    d_tocheck = d_record(tocheck,:);
    [f_val_aux, u_val_record_aux, nfeval_aux] = u_validation(problem_max_u, d_tocheck, u_record, lsflag_validation, 1:n_obj);
    nfeval_val = nfeval_val + nfeval_aux;

    % update f_record and u_val_record
    f_record(tocheck,:) = f_val_aux;
    for obj = 1:n_obj
        u_val_record{obj}(tocheck,:) = u_val_record_aux{obj};
    end
    checked = checked | tocheck;           % update who has been checked
end

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

