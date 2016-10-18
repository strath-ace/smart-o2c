function [d,fval,exitflag,output] = optimise_minmarek(problem_minmax, algo_outer, algo_inner, par_minmax)

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
tol_conv = par_minmax.tol_conv; % Tolerance to convergence, only for S.O.

lsflag_validation = par_minmax.local_search_flags.validation; % run either local search or func eval in validation
keep_d_record = par_minmax.keep_d_record;

% Hardcoded MinMaReK inputs
n_d_0 = 10*n_d; %HC
n_u_0 = 10*n_u; %HC
n_u_record_0 = repmat([1],[1,n_obj]); %HC

% Build a mock macsminmax metaproblem max_u that we will use for validation and function evaluation
problem_max_u = build_metaproblem_macsminmax_inner(problem_minmax);
objfun = problem_max_u.objfun;

% Build true inner and outer problems of MinMaReK
% This will also initialise the surrogates
problem_outer = build_metaproblem_minmarek_outer(problem_minmax);
problem_inner = build_metaproblem_minmarek_inner(problem_minmax);

% Random initial guesses for the DOE in d and u and for the archive in u
d_0 = lhsu(zeros(1,n_d),ones(1,n_d),n_d_0);
u_0 = lhsu(zeros(1,n_u),ones(1,n_u),n_u_0);
u_record = cell(1,n_obj);
for obj = 1:n_obj
    u_record{obj} = lhsu(zeros(1,n_u),ones(1,n_u),n_u_record_0(obj));
end
d_record = [];
f_record = [];

% related to saving fevals thanks to the fact that we recycle the DOE in d (d_0)
u_record_aux = u_record;
f_d_0 = -sign_inner*inf(n_d_0,n_obj);


% Now the fun begins
nfeval = 0;

% Main loop
iter=0;
precision = inf;
stop = false;
while ~stop

    iter = iter +1;

    %OUTER LOOP: MINIMISATION OVER D
    % related to saving fevals by recycling d_0
    [f_d_0_aux, ~, nfeval_aux] = u_validation(problem_max_u, d_0, u_record_aux, false, 1:n_obj); %here we can save some fevals with a val_record
    nfeval = nfeval + nfeval_aux;

    u_record_aux = cell(1,n_obj);
    f_d_0 = sign_inner*max(sign_inner*f_d_0, sign_inner*f_d_0_aux);

    f_outer = f_d_0;
    d_outer = d_0;
    
    indicator_d = inf;
    iter_d = 0;
    while abs(indicator_d) >= indicator_d_max && iter_d < iter_d_max

        iter_d = iter_d+1;

        % update surrogate and ymin
        [~,idx] = unique(round(1e8*d_outer),'rows');
        d_outer = d_outer(idx,:);
        f_outer = f_outer(idx,:);

        problem_outer.par_objfun.surrogate.model = problem_outer.par_objfun.surrogate.training(d_outer,f_outer,problem_outer.par_objfun.surrogate);
        if n_obj == 1
            problem_outer.par_objfun.ymin = min(f_outer);
        else
            problem_outer.par_objfun.ymin = f_outer(dominance(f_outer,0) == 0,:); % pareto front
        end

        %maximise indicator_d
        [d_outer_aux, indicator_d, ~, ~] = algo_outer.optimise(problem_outer,algo_outer.par); %no feval here
        indicator_d = -indicator_d;

        %validate the result
        [f_outer_aux,~,nfeval_aux] = u_validation(problem_max_u, d_outer_aux, u_record, false, 1:n_obj);
        nfeval = nfeval + nfeval_aux;

        %apend d and f
        d_outer = [d_outer; d_outer_aux];
        f_outer = [f_outer; f_outer_aux];
        
    end

    iter_d
    indicator_d

    %% Select non-dominated solutions
    sel = dominance(f_outer,0) == 0;
    fmin = f_outer(sel,:);
    dmin = d_outer(sel,:);
    
    % I don't think this is necessary but...
    dmin(dmin < 0.0) = 0.0;
    dmin(dmin > 1.0) = 1.0;

    % Remove solutions archived more than once (if any)
    [~,idx] = unique(round(1e8*dmin),'rows');
    dmin = dmin(idx,:);
    fmin = fmin(idx,:);
    n_dmin = size(dmin,1); 
    
    if (n_obj == 2)
        figure(1)
        colors = 'ckbrmg';
        hold on
        plot(fmin(:,1),fmin(:,2),strcat(colors(mod(iter,length(colors))+1),'.'))
        figure(1)    
    end
    fmin

    % Archive shrinking for MO?
    
    %Convergence
    if n_obj == 1
        precision = fmin;
    end
    
    % INNER LOOP: MAXIMISATION OVER U
    fmax_inner = fmin;
    for i = 1:n_dmin
        problem_max_u.par_objfun.d = dmin(i,:);

        for obj = 1:n_obj
            problem_max_u.par_objfun.objective = obj;

            u_inner = u_0;
            f_inner =[];
            for j = 1:size(u_inner,1);
                f_inner(j,1)=objfun(u_inner(j,:),problem_max_u.par_objfun);
                nfeval = nfeval + 1;
            end
            
            indicator_u = inf;
            iter_u = 0;
            while abs(indicator_u) >= indicator_u_max && iter_u < iter_u_max
                iter_u = iter_u + 1;

                % update surrogate and ymin
                [~,idx] = unique(round(1e8*u_inner),'rows');
                u_inner = u_inner(idx,:);
                f_inner = f_inner(idx,:);
                problem_inner.par_objfun.surrogate.model = problem_inner.par_objfun.surrogate.training(u_inner,f_inner,problem_inner.par_objfun.surrogate);
                problem_inner.par_objfun.ymin = min(f_inner);

                % maximise indicator_u
                [u_inner_aux, indicator_u, ~, ~] = algo_inner.optimise(problem_inner,algo_inner.par); %no nfeval here
                indicator_u = -indicator_u;

                f_inner_aux = objfun (u_inner_aux, problem_max_u.par_objfun); %evaluate. note it returns negative for minmax
                nfeval = nfeval + 1;
                
                %apend u and f
                u_inner = [u_inner; u_inner_aux];
                f_inner = [f_inner; f_inner_aux];
                
            end

            iter_u
            indicator_u

            [fmax, sel] = min(f_inner);
            fmax = -sign_inner*fmax;
            u_record_aux{obj} = [u_record_aux{obj}; u_inner(sel,:)];

            if sign_inner*fmax > sign_inner*fmax_inner(i,obj)
                fmax_inner(i,obj) = fmax;
            end
        end
        
        if (n_obj == 2)
            figure(1)
            colors = 'ckbrmg';
            hold on
            plot(fmax_inner(i,1),fmax_inner(i,2),strcat(colors(mod(iter,length(colors))+1),'o'))
            figure(1)
        end
        [fmin(i,:) ; fmax_inner(i,:)]
    end

    if (n_obj == 2)
        figure(1)
        colors = 'ckbrmg';
        hold on
        plot(fmax_inner(:,1),fmax_inner(:,2),strcat(colors(mod(iter,length(colors))+1),'o'))
        figure(1)
    end
    fmax_inner

    d_record = [d_record;dmin];
    f_record = [f_record;fmax_inner];

    for obj = 1:n_obj
        u_record{obj} = [u_record{obj}; u_record_aux{obj}];
        [~,idx] = unique(round(1e8*u_record{obj}),'rows');
        u_record{obj} = u_record{obj}(idx,:);
    end

    size_u_record = sum(cellfun('size',u_record,1));
    nfeval_val = size_u_record*size(dmin,1);% (keep_d_record*size(d_record,1) + ~keep_d_record*size(dmin,1)); % Quite imprecise, should add rule of thumb
    stop = nfeval + nfeval_val >= nfevalmax;
    if n_obj == 1
        precision = sign_inner*(fmax - precision);
        stop = stop || precision < tol_conv;
    end
    
    nfeval
    nfevalglobal
end

if ~keep_d_record
    d_record = dmin;
    f_record = fmax_inner;
end

[~,idx] = unique(round(1e8*d_record),'rows');
d_record = d_record(idx,:);
f_record = f_record(idx,:);

%% more clever archive cross check
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
exitflag = nfeval <= nfevalmax && precision < tol_conv; 

end

