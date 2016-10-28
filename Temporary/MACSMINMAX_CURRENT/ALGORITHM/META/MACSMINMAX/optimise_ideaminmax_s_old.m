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
x_doe = cell(1,n_obj);
f_doe = cell(1,n_obj);
for i=1:n_du_0
    problem_max_u.par_objfun.d = du_0(i,1:n_d);
    for obj = 1:n_obj
        problem_max_u.par_objfun.objective = obj;

        x_doe{obj}(i,:) = du_0(i,:);
        f_doe{obj}(i,:) = -sign_inner*problem_max_u.objfun(x_doe{obj}(i,n_d+1:n_d+n_u),problem_max_u.par_objfun);
    end
    nfeval = nfeval + 1;
end

for obj = 1:n_obj
%     problem_outer.par_objfun.surrogate.model{obj} = problem_outer.par_objfun.surrogate.training(x_doe{obj},f_doe{obj},problem_outer.par_objfun.surrogate);
    problem_inner.par_objfun.surrogate.model{obj} = problem_inner.par_objfun.surrogate.training(x_doe{obj},-sign_inner*f_doe{obj},problem_inner.par_objfun.surrogate);
end


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
    
    [~, idmaxdoe] = max(f_doe{obj});
    dmaxdoe = x_doe{obj}(idmaxdoe,1:n_d);
    
    problem_max_u.par_objfun.d = dmaxdoe;
    problem_max_u.par_objfun.objective = obj;
    
    problem_inner.par_objfun.d = dmaxdoe;
    problem_inner.par_objfun.surrogate.model{obj} = problem_inner.par_objfun.surrogate.training(x_doe{obj},-sign_inner*f_doe{obj},problem_outer.par_objfun.surrogate);
    
    [umaxmaxdoe , ~ , ~ , ~ ] = algo_inner.optimise(problem_inner,algo_inner.par);
    f_aux = -sign_inner*problem_max_u.objfun(umaxmaxdoe,problem_max_u.par_objfun);
    nfeval = nfeval + 1;
    
    u_record{obj}=[u_record{obj};umaxmaxdoe];
    
    x_doe{obj} = [x_doe{obj};[dmaxdoe umaxmaxdoe]];
    f_doe{obj} = [f_doe{obj};f_aux];
    
    problem_outer.par_objfun.surrogate.model{obj} = problem_outer.par_objfun.surrogate.training(x_doe{obj},f_doe{obj},problem_outer.par_objfun.surrogate);
end


% Main loop
iter=0;
precision = inf;
problem_outer.par_objfun.ymin = 1e32;
problem_inner.par_objfun.ymin = 1e32;
dmin_outer = nan(1,n_d);
fmin_outer = 1e32;
umax_inner = nan(1,n_d);
fmax_inner = sign_inner*1e32;
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
        % maximise indicator_d
        [d_outer_aux, indicator_d, ~, ~] = algo_outer.optimise(problem_outer,algo_outer.par); %no feval here
        indicator_d = -indicator_d;

        % validate dmin with the archive
        [f_outer_aux, u_outer_aux, nfeval_aux, all_f] = u_validation(problem_max_u, d_outer_aux, u_record, false, 1:n_obj);
        nfeval = nfeval + nfeval_aux;

        % update x_doe, f_doe and surrogate
        % NOTE suppose there is one doe per objective and we'll see if this works or not...
        % NOTE II in fact this is stupid because evaluating one objective should be more or less as hard as evaluating them all... 
        % so I should recode with just one doe asap ("recoding" validation here) but for the moment leave it like this for SO version
        for obj = 1:n_obj
            size_u_record_obj = size(u_record{obj},1);
            for i = 1:size(d_outer_aux,1)
                x_doe{obj} = [x_doe{obj}; repmat(d_outer_aux(i,:),[size_u_record_obj,1]), u_record{obj}];
                f_doe{obj} = [f_doe{obj}; all_f{obj}(i,:)'];
            end
            % [x_doe{obj},idx] = uniquetol(x_doe{obj},1e-8,'ByRows',true); % R2015+
            [~, idx] = unique(round(1e8*x_doe{obj}),'rows');
            x_doe{obj} = x_doe{obj}(idx,:);
            f_doe{obj} = f_doe{obj}(idx,:);
            problem_outer.par_objfun.surrogate.model{obj} = problem_outer.par_objfun.surrogate.training(x_doe{obj},f_doe{obj},problem_outer.par_objfun.surrogate);
        end

        % update fmin_outer and friends, and ymin for EI computation
        if n_obj == 1
            [fmin_outer,idx] = min([fmin_outer;f_outer_aux]);
            dmin_outer = [dmin_outer; d_outer_aux];
            dmin_outer = dmin_outer(idx,:);
            problem_outer.par_objfun.ymin = fmin_outer; % we'll come back to this, but it's like this in ideaminmax
        else
            error('MO: some things not implemented! (and nothing tested...)')
            % problem_outer.par_objfun.ymin = f_outer(dominance(f_outer,0) == 0,:); % pareto front
            % probably the simpler is to make a stack f_outer with all f_outer_aux and choose the non-dominated
        end
        
    end

% %% TAKE CARE LATER OF ALL OF THIS (UNTIL NEXT WHILE LOOP)
% %% REMEMBER TO TRAIN THE SURROGATE FOR THE INNER LOOP HERE
%     %% Select non-dominated solutions
%     sel = dominance(f_outer,0) == 0;
%     fmin = f_outer(sel,:);
%     dmin = d_outer(sel,:);
    
%     % I don't think this is necessary but...
%     dmin(dmin < 0.0) = 0.0;
%     dmin(dmin > 1.0) = 1.0;

%     % Remove solutions archived more than once (if any)
%     [dmin,idx] = unique(dmin,'rows');
%     fmin = fmin(idx,:);
    n_dmin = size(dmin_outer,1);
    % problem_inner.par_objfun.surrogate.model = problem_outer.par_objfun.surrogate.model;
    for obj = 1:n_obj
        problem_inner.par_objfun.surrogate.model{obj} = problem_inner.par_objfun.surrogate.training(x_doe{obj},-sign_inner*f_doe{obj},problem_inner.par_objfun.surrogate);
    end
    % % Archive shrinking for MO??
    
    % INNER LOOP: MAXIMISATION OVER U
    for i = 1:n_dmin
        problem_inner.par_objfun.d = dmin_outer(i,:);
        problem_max_u.par_objfun.d = dmin_outer(i,:);
        f_record_aux=zeros(1,n_obj);
        for obj = 1:n_obj
            problem_max_u.par_objfun.objective = obj;
            
            indicator_u = inf;
            iter_u = 0;
            while abs(indicator_u) > indicator_u_max && iter_u < iter_u_max
                iter_u = iter_u + 1;

                % maximise indicator_u
                [u_inner_aux, indicator_u, ~, ~] = algo_inner.optimise(problem_inner,algo_inner.par); %no nfeval here
                indicator_u = -indicator_u;

                % evaluate
                f_inner_aux = problem_max_u.objfun(u_inner_aux,problem_max_u.par_objfun); %note it comes negative for minmax
                nfeval = nfeval + 1;

                % update doe and surrogate
                % if ~ismembertol([dmin_outer(i,:) , u_inner_aux],x_doe{obj},1e-8,'ByRows',true) %R2015+
                if ~ismember(round(1e8*[dmin_outer(i,:) , u_inner_aux]),round(1e8*x_doe{obj}),'rows') 
                    x_doe{obj} = [x_doe{obj}; dmin_outer(i,:) , u_inner_aux];
                    f_doe{obj} = [f_doe{obj}; -sign_inner*f_inner_aux];
                end %should be enough to substitute the more expensive 'unique'
                problem_inner.par_objfun.surrogate.model{obj} = problem_inner.par_objfun.surrogate.training(x_doe{obj},-sign_inner*f_doe{obj},problem_inner.par_objfun.surrogate);

                % update ymin for EI computation
                if n_obj == 1
                    [fmax_inner,idx] = min([fmax_inner;f_inner_aux]);
                    umax_inner = [umax_inner; u_inner_aux];
                    umax_inner = umax_inner(idx,:);
                    problem_inner.par_objfun.ymin = fmax_inner; % we'll come back to this, but it's like this in ideaminmax
                else
                    error('MO: some things not implemented! (and nothing tested...)')
                    % problem_outer.par_objfun.ymin = f_outer(dominance(f_outer,0) == 0,:); % pareto front
                end
            end

            if fmax_inner < -sign_inner*fmin_outer(i,obj)
                u_record_aux{obj} = [u_record_aux{obj}; umax_inner];
                f_record_aux(1,obj) = -sign_inner*fmax_inner;
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
    
    
    % d_record = uniquetol(d_record,1e-8,'ByRows',true); % R2015+
    [~,idx] = unique(round(1e8*d_record),'rows');
    d_record = d_record(idx,:);
    f_record = f_record(idx,:); % not worrying too much if i lost validated values
    for obj = 1:n_obj
        if iter==1
            u_record{obj} = []; % I don't like this too much
        end
        u_record{obj} = [u_record{obj}; u_record_aux{obj}];
        % u_record{obj} = uniquetol(u_record{obj},1e-8,'ByRows',true); % R2015b+
        [~,idx] = unique(round(1e8*u_record{obj}),'rows');
        u_record{obj} = u_record{obj}(idx,:);

        problem_outer.par_objfun.surrogate.model{obj} = problem_outer.par_objfun.surrogate.training(x_doe{obj},f_doe{obj},problem_outer.par_objfun.surrogate);
    end
    
    size_u_record = sum(cellfun('size',u_record,1));
    nfeval_val = size_u_record*size(dmin_outer,1);
    stop = nfeval + nfeval_val >= nfevalmax;
    if n_obj == 1
        % precision = sign_inner*(-sign_inner*fmax_inner - fmin_outer);
        precision = abs(-sign_inner*fmax_inner - fmin_outer);
        stop = stop || precision < tol_conv;
    end
    
    fmax_inner = -sign_inner*min(sign_inner*f_doe{1});
    fmin_outer = max(f_doe{1});
    
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
nfeval_val = 0;
while ~stop
    sel = dominance(f_record,0) == 0;       % find non-dominated
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

