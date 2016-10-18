function [algo_minmax, algo_outer, algo_inner] = init_ideaminmax_s(minmax_problem, varargin);

%% META ALGORITHM: MACSMINMAX
algo_minmax.optimise = @optimise_ideaminmax_s;            % algorithm in the form [x,f,exitflag,output] = algo(problem,algo_outer,algo_inner,par_minmax)
% parameters
    if isfield(minmax_problem,'maxnfeval')
        par_minmax.maxnfeval = minmax_problem.maxnfeval/10;% function evaluation limit (does not account for ls in validation)
    else
        error('max nfeval not supplied')
    end

    par_minmax.local_search_flags.validation = true;
    par_minmax.keep_d_record = true;
    par_minmax.tol_conv = 1e-3; %Tolerance to convergence, only used in S.O. case
    par_minmax.indicator_d_max = 1e-3; % Threshold on the indicator (EI, PI, etc.) for outer loop
    par_minmax.indicator_u_max = 1e-3; % Threshold on the indicator (EI, PI, etc.) for inner loop
    par_minmax.iter_d_max = 20; %20*minmax_problem.dim_d; % Max. iterations for outer loop in an iteration of MinMaReK
    par_minmax.iter_u_max = 20; %20*minmax_problem.dim_u; % Max. iterations for inner loop in an iteration of MinMaReK


algo_minmax.par_minmax = par_minmax;


%% ALGORITHM OUTER LOOP: IDEA2
algo_outer.optimise = @optimise_idea2_wrapper;            % algorithm in the form [x,f,exitflag,output] = algo(problem,par_algo)
    %parameters{
        options_idea2(1) = 150*minmax_problem.dim_d;                             % number of function evaluations for IDEA
        options_idea2(4) = 5;                                 % number of agents
        options_idea2(31) = 1;                                % number of local restarts
        options_idea2(32) = 1;                                % communication heuristics (see com_red.m)
        options_idea2(22) = 0.1;                              % delta_c: crowding factor for global restart
        options_idea2(27) = 0.25;                             % tolconv: local convergence
        options_idea2(15) = -1;                               % verbousity
        options_idea2(28) = 0.2;                                % F
        options_idea2(29) = 0.8;                                % CR
        % options_idea2(30) = minmax_problem.n_obj;             % number of objectives (not used)
        options_idea2(33) = 0.2;                              % initial probability of running the subproblem (not used)

algo_outer.par.options = options_idea2;
algo_outer.par.initial_population = []; 


%% ALGORITHM INNER LOOP: IDEA2
algo_inner.optimise = @optimise_idea2_wrapper;            % algorithm in the form [x,f,exitflag,output] = algo(problem,par_algo)
    %parameters{
        options_idea2(1) = 150*minmax_problem.dim_u;          % number of function evaluations for IDEA
        options_idea2(4) = 5;                                 % number of agents
        options_idea2(31) = 1;                                % number of local restarts
        options_idea2(32) = 1;                                % communication heuristics (see com_red.m)
        options_idea2(22) = 0.1;                              % delta_c: crowding factor for global restart
        options_idea2(27) = 0.25;                             % tolconv: local convergence
        options_idea2(15) = -1;                               % verbousity
        options_idea2(28) = 0.2;                                % F
        options_idea2(29) = 0.8;                                % CR
        % options_idea2(30) = minmax_problem.n_obj;             % number of objectives (not used)
        options_idea2(33) = 0.2;                              % initial probability of running the subproblem (not used)

algo_inner.par.options = options_idea2;
algo_inner.par.initial_population = [];

return
