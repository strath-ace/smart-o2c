function [algo_inner] = init_algorithm(problem)
%% ALGORITHM INNER LOOP: MPAIDEA
algo_inner.optimise = @optimise_mpaidea_wrapper_decomposition;      % algorithm in the form [x,f,exitflag,output] = algo(problem,par_algo)
par_mpaidea.nFeValMax = 100;%*problem.dim_u;            % number of function evaluations for IDEA
par_mpaidea.n_populations = 1;                        % number of populations, if no adaptive behaviour should set to 1
par_mpaidea.n_agents = max(5,problem.dim_u);          % number of agents in one population
par_mpaidea.population=[];                            % initial population, same for all execution. Leave empty to randomize.
par_mpaidea.max_LR = 100;                             % max.number of local restarts
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
end