function [algo_minmax, algo_outer, algo_inner] = init_algo_so(minmax_problem)

%% META ALGORITHM: MACSMINMAX
algo_minmax.optimise = @optimise_so;            % algorithm in the form [x,f,exitflag,output] = algo(problem,algo_outer,algo_inner,par_minmax)
% parameters
    if isfield(minmax_problem,'maxnfeval')
        par_minmax.maxnfeval = minmax_problem.maxnfeval;% TOTAL function evaluation limit
    else
        par_minmax.maxnfeval = 1e4;                     %error('max nfeval not supplied') % 100;%
    end
    par_minmax.n_d0 = 1;
    
    par_minmax.local_search_flags.validation = false;    % Run local search (true) or function evaluation (false) in validation
    par_minmax.local_search_flags.inner = false;        % Idem inner (SO) problem
    par_minmax.local_search_flags.outer = false;       % Idem outer (MO) problem
    % par_minmax.tconv = 1;                               % tconv
    % par_minmax.tspr = 1;                                % tspr
algo_minmax.par_minmax = par_minmax;


if (minmax_problem.n_obj == 1)
    %% ALGORITHM OUTER LOOP: MPAIDEA
    algo_outer.optimise = @optimise_mpaidea_wrapper;          % algorithm in the form [x,f,exitflag,output] = algo(problem,par_algo)
        par_mpaidea.nFeValMax = 5000;%1000*minmax_problem.dim_d;     % number of function evaluations for IDEA     100;%
                    par_mpaidea.n_populations = 4;                        % number of populations, if no adaptive behaviour should set to 1
        par_mpaidea.n_agents = max(5,minmax_problem.dim_d);   % number of agents in one population
        par_mpaidea.population=[];                            % initial population, same for all execution. Leave empty to randomize.
                    par_mpaidea.max_LR = [];                             % max.number of local restarts
        par_mpaidea.DE_strategy = 1;                          % 1: DE/Rand and DE/CurrentToBest; 2: DE/Rand and DE/Best
        par_mpaidea.prob_DE_strategy = 0.5;                   % probability or not using DE/Rand
        par_mpaidea.delta_local = 0.1;                        % dimension of the bubble for the local restart of population.
        par_mpaidea.delta_global = 0.1;                       % characteristic dimension for the global restart of the population
        par_mpaidea.rho = 0.25;                               % contraction threshold for the population
        par_mpaidea.dd_CRF = 3;                               % parameter for the adaptation of CRF.
        par_mpaidea.F = 0.2;                                  % F
        par_mpaidea.CR = 0.8;                                 % CR
                    par_mpaidea.text = 0;                                 % verbosity
        par_mpaidea.record = [1];                             % Fraction(s) of nfevalmax at which the result is saved
    algo_outer.par = par_mpaidea;


    % %% ALGORITHM OUTER LOOP: IDEA2
    % algo_outer.optimise = @optimise_idea2_wrapper;            % algorithm in the form [x,f,exitflag,output] = algo(problem,par_algo)
    % %parameters{
    %     options_idea2(1) = 200*minmax_problem.dim_u;          % number of function evaluations for IDEA
    %     options_idea2(4) = 5;                                 % number of agents
    %     options_idea2(31) = 100;                              % number of local restarts
    %     options_idea2(32) = 2;                                % communication heuristics (see com_red.m)
    %     options_idea2(22) = 0.1;                              % delta_c: crowding factor for global restart
    %     options_idea2(27) = 0.25;                             % tolconv: local convergence
    %     options_idea2(15) = -1;                               % verbousity
    %     options_idea2(28) = 0.2;                                % F
    %     options_idea2(29) = 0.8;                                % CR
    %     % options_idea2(30) = minmax_problem.n_obj;             % number of objectives (not used)
    %     options_idea2(33) = 0.2;                              % initial probability of running the subproblem (not used)

    % algo_outer.par.options = options_idea2;
    % algo_outer.par.initial_population = []; 
else
    warning('using SO settings for a MO problem!')
    %% ALGORITHM OUTER LOOP: MACS
    algo_outer.optimise = @optimise_macs_wrapper;           % algorithm in the form [x,f,exitflag,output] = algo(problem,par_algo)
    %parameters{
        par_macs.maxnfeval = minmax_problem.dim_d;      % Max function evaluations   100;%
        par_macs.popsize = 1;                              % Population size
        par_macs.rhoini = 1;                                % Initial local hypercube size
        par_macs.F = 1;                                     % F
        par_macs.CR = 0.1;                                  % CR
        par_macs.p_social = 0.5;                            % Ratio between elite and total population
        par_macs.max_arch = 20;                             % Output archive size
        par_macs.coord_ratio = 1;                           % Quota of coordinates examined for each set of individualistic actions
        par_macs.contr_ratio = 0.5;                         % Contraction ratio
        par_macs.draw_flag = 0;                             % Print itarations status  
        par_macs.cp = 0;                                    % constraint flag
        par_macs.MBHflag = 0;                               % number of MBH steps
    algo_outer.par = par_macs;
end

%% ALGORITHM INNER LOOP: MPAIDEA
algo_inner.optimise = @optimise_mpaidea_wrapper;          % algorithm in the form [x,f,exitflag,output] = algo(problem,par_algo)
    par_mpaidea.nFeValMax = 5000;%100*minmax_problem.dim_u;     % number of function evaluations for IDEA    100;%
                    par_mpaidea.n_populations = 4;                        % number of populations, if no adaptive behaviour should set to 1
    par_mpaidea.n_agents = max(5,minmax_problem.dim_u);   % number of agents in one population
    par_mpaidea.population=[];                            % initial population, same for all execution. Leave empty to randomize.
                    par_mpaidea.max_LR = [];                             % max.number of local restarts
    par_mpaidea.DE_strategy = 1;                          % 1: DE/Rand and DE/CurrentToBest; 2: DE/Rand and DE/Best
    par_mpaidea.prob_DE_strategy = 0.5;                   % probability or not using DE/Rand
    par_mpaidea.delta_local = 0.1;                        % dimension of the bubble for the local restart of population.
    par_mpaidea.delta_global = 0.1;                       % characteristic dimension for the global restart of the population
    par_mpaidea.rho = 0.25;                               % contraction threshold for the population
    par_mpaidea.dd_CRF = 3;                               % parameter for the adaptation of CRF.
    par_mpaidea.F = 0.2;                                  % F
    par_mpaidea.CR = 0.8;                                 % CR
    par_mpaidea.text = 0;                                 % verbosity
    par_mpaidea.record = [1];                             % Fraction(s) of nfevalmax at which the result is saved
algo_inner.par = par_mpaidea;




% %% ALGORITHM INNER LOOP: IDEA2
% algo_inner.optimise = @optimise_idea2_wrapper;            % algorithm in the form [x,f,exitflag,output] = algo(problem,par_algo)
% %parameters{
%     options_idea2(1) = 200*minmax_problem.dim_u;          % number of function evaluations for IDEA
%     options_idea2(4) = 5;                                 % number of agents
%     options_idea2(31) = 100;                              % number of local restarts
%     options_idea2(32) = 2;                                % communication heuristics (see com_red.m)
%     options_idea2(22) = 0.1;                              % delta_c: crowding factor for global restart
%     options_idea2(27) = 0.25;                             % tolconv: local convergence
%     options_idea2(15) = -1;                               % verbousity
%     options_idea2(28) = 0.2;                                % F
%     options_idea2(29) = 0.8;                                % CR
%     % options_idea2(30) = minmax_problem.n_obj;             % number of objectives (not used)
%     options_idea2(33) = 0.2;                              % initial probability of running the subproblem (not used)

% algo_inner.par.options = options_idea2;
% algo_inner.par.initial_population = [];
return