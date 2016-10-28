function [algo_minmax, algo_outer, algo_inner] = init_mo_tc2(minmax_problem);

%% META ALGORITHM: MACSMINMAX
algo_minmax.optimise = @optimise_mo;            % algorithm in the form [x,f,exitflag,output] = algo(problem,algo_outer,algo_inner,par_minmax)
% parameters
    if isfield(minmax_problem,'maxnfeval')
        par_minmax.maxnfeval = minmax_problem.maxnfeval*20;% TOTAL function evaluation limit
    else
        error('max nfeval not supplied')
    end
    par_minmax.n_d0 = 8;

    par_minmax.local_search_flags.validation = true;    % Run local search (true) or function evaluation (false) in validation
    par_minmax.local_search_flags.inner = false;        % Idem inner (SO) problem
    par_minmax.local_search_flags.outer = false;        % Idem outer (MO) problem
    % par_minmax.tconv = 1;                               % tconv
    % par_minmax.tspr = 1;                                % tspr
algo_minmax.par_minmax = par_minmax;

if (minmax_problem.n_obj == 1)
    %% ALGORITHM OUTER LOOP: IDEA2
    algo_outer.optimise = @optimise_idea2_wrapper;            % algorithm in the form [x,f,exitflag,output] = algo(problem,par_algo)
    %parameters{
        options_idea2(1) = 200*minmax_problem.dim_u;          % number of function evaluations for IDEA
        options_idea2(4) = 5;                                 % number of agents
        options_idea2(31) = 20;                               % number of local restarts
        options_idea2(32) = 2;                                % communication heuristics (see com_red.m)
        options_idea2(22) = 0.1;                              % delta_c: crowding factor for global restart
        options_idea2(27) = 0.1;                              % tolconv: local convergence
        options_idea2(15) = -1;                               % verbousity
        options_idea2(28) = 1;                                % F
        options_idea2(29) = 0.1;                              % CR
        % options_idea2(30) = minmax_problem.n_obj;             % number of objectives (not used)
        options_idea2(33) = 0.2;                              % initial probability of running the subproblem (not used)

    algo_outer.par.options = options_idea2;
    algo_outer.par.initial_population = []; 
else
    %% ALGORITHM OUTER LOOP: MACS
    algo_outer.optimise = @optimise_macs_wrapper;           % algorithm in the form [x,f,exitflag,output] = algo(problem,par_algo)
    %parameters{
        par_macs.maxnfeval = 400*minmax_problem.dim_d;      % Max function evaluations
        par_macs.popsize = 3*minmax_problem.dim_d;          % Population size
        par_macs.rhoini = 1;                                % Initial local hypercube size
        par_macs.F = 0.9;                                   % F
        par_macs.CR = 0.9;                                  % CR
        par_macs.p_social = 0.2;                            % Ratio between elite and total population
        par_macs.max_arch = 40;                             % Internal archive size in MACS
        par_macs.max_arch_out = 8;                          % Output archive size.
        par_macs.coord_ratio = 1;                           % Quota of coordinates examined for each set of individualistic actions
        par_macs.contr_ratio = 0.5;                         % Contraction ratio
        par_macs.draw_flag = 0;                             % Print itarations status  
        par_macs.cp = 0;                                    % constraint flag
        par_macs.MBHflag = 0;                               % number of MBH steps
    algo_outer.par = par_macs;
end

%% ALGORITHM INNER LOOP: IDEA2
algo_inner.optimise = @optimise_idea2_wrapper;            % algorithm in the form [x,f,exitflag,output] = algo(problem,par_algo)
%parameters{
    options_idea2(1) = 400*minmax_problem.dim_u;          % number of function evaluations for IDEA
    % options_idea2(1) = 500;
    options_idea2(4) = 2*minmax_problem.dim_u;            % number of agents
    options_idea2(31) = 20;                               % number of local restarts
    options_idea2(32) = 2;                                % communication heuristics (see com_red.m)
    options_idea2(22) = 0.1;                              % delta_c: crowding factor for global restart
    options_idea2(27) = 0.1;                              % tolconv: local convergence
    options_idea2(15) = -1;                               % verbousity
    options_idea2(28) = .7;                               % F
    options_idea2(29) = .6  ;                             % CR
    % options_idea2(30) = minmax_problem.n_obj;           % number of objectives (not used)
    options_idea2(33) = 0.2;                              % initial probability of running the subproblem (not used)

algo_inner.par.options = options_idea2;
algo_inner.par.initial_population = [];

return