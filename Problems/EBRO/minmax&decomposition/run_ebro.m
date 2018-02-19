% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example of run of ebro problem 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
close all
clc


global nfevalglobal;
nfevalglobal = 0;


% Reset random numbers generator
seed = 1;
s = RandStream('mt19937ar','Seed', seed);
RandStream.setGlobalStream(s);


% Add path to optimiser folder
addpath(genpath('Optimisation'))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Choose the output

% -------------------------------------------------------------------------
% problem.output = 0 --> only minmax-Belief
% problem.output = 1 --> only minmin-Plausibility
% problem.output = 2 --> both minmax-Belief and minmin-Plaiusibility
% -------------------------------------------------------------------------
problem.output = 0;


% -------------------------------------------------------------------------
% problem.input = 0 --> do minmax and minmin
% problem.input = 1 --> load d, u_min, u_max
% problem.input = 2 --> load d, run max and min
% -------------------------------------------------------------------------
problem.input = 0;


% -------------------------------------------------------------------------
% if input == 1 or input == 2.
% if you have the structure(s) 'minmax' or 'minmin' from the 
% optimise_macs_minmax problem write the name of the file .mat here:
% -------------------------------------------------------------------------
problem.input_minmin_minmax = 'minmax_minmin_TC1.mat'; 

% -------------------------------------------------------------------------
% oterwise problem.input.problem.input.minmin_minmax = [] and:
% -------------------------------------------------------------------------
problem.input_d_min = [];    % d_min    = [d1, d2, ..., dn] 
problem.input_d_max = [];    % d_max    = [d1, d2, ..., dn] 
problem.input_u_min = [];    % u_min{1} = [u1, u2, ...,un]
problem.input_u_max = [];    % u_max{1} = [u1, u2, ...,un]
problem.input_F_min = [];    % F_min    = F 
problem.input_F_max = [];    % F_max    = F





% -------------------------------------------------------------------------
% do decompoition and/or exact Belief
% problem.exact_curves = 0 --> no
% problem.exact_curves = 1 --> yes
% problem.exact_curves = 2 --> DO ONLY EXACT CURVE
% -------------------------------------------------------------------------
problem.exact_curves = 1;


% number of sub-functions decomposition
num_functions = 2;  % number of sub-functions in which the problem is decomposed


% number of samples
num_samples = 5;    % number of samples for each Belief and Plausibility curve of coupled vector


problem.num_functions = num_functions;
for i = 1:problem.num_functions/2*(problem.num_functions-1)
    problem.num_samples{i} = num_samples;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Problem input

%--------------------------------------------------------------------------
% Dimention, Lower and Upper boundaries of the EPISTEMIC parameters 
%--------------------------------------------------------------------------
problem.dim_u = [2 2 2]; % [u1, u2, u12]

% 7) bounds of uncertain vector u;
problem.lb_u{1} = {[-5 -3 0]; [-5 -3 0];...
    [-5 -3 0]; [-5 -3 0];...
    [-5 -3 0]; [-5 -3 0]};

problem.ub_u{1} = {[-1 0 2]; [-1 0 2];...
    [-1 0 2]; [-1 0 2];...
    [-1 0 2]; [-1 0 2]};

problem.bpa{1} =  {[.3 .3 .4]; [.3 .3 .4];...
    [.3 .3 .4]; [.3 .3 .4];...
    [.3 .3 .4]; [.3 .3 .4]};

%--------------------------------------------------------------------------
% Dimention, Lower and Upper boundaries of the DESIGN parameters 
%--------------------------------------------------------------------------
problem.dim_d = 2;
problem.lb_d = [-5; -5];
problem.ub_d = [5; 5];


%--------------------------------------------------------------------------
% List of the FIXED parameters
%--------------------------------------------------------------------------

problem.fix = [];

problem.fix = problem.fix;
problem.par_objfun{1}.fix = problem.fix;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Objective function


%--------------------------------------------------------------------------
% Number of objective functions
%--------------------------------------------------------------------------
problem.n_obj = 1;

%--------------------------------------------------------------------------
% Objective functions;
%--------------------------------------------------------------------------
problem.objfun =     {@TC_1};

%--------------------------------------------------------------------------
% Constraints
%--------------------------------------------------------------------------
% problem.constraints = {[]}; 
problem.constraints = {@TC_1_constraints};



% design variables
problem.dim_d = problem.dim_d;
problem.lb_d = problem.lb_d;
problem.ub_d = problem.ub_d;
% uncertain variables
problem.dim_u_i = problem.dim_u;
dim_u = sum(problem.dim_u_i);
problem.dim_u = dim_u;

for n =1:problem.n_obj
    problem.lb_u{n} = problem.lb_u{n};
    problem.ub_u{n} = problem.ub_u{n};
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MP_AIDEA Parameters


% META LOOP

%--------------------------------------------------------------------------
% algorithm in the form [x,f,exitflag,output] = algo(problem,algo_outer,algo_inner,par_minmax)
%--------------------------------------------------------------------------
algo_minmax.optimise = @optimise_so;  

%--------------------------------------------------------------------------
% maximum number of function evaluation  
%--------------------------------------------------------------------------
    par_minmax.maxnfeval = 1e5;   

%--------------------------------------------------------------------------
% number of initial design vectors
%--------------------------------------------------------------------------
    par_minmax.n_d0 = 1;

%--------------------------------------------------------------------------
% Run local search (true) or function evaluation (false) in validation
%--------------------------------------------------------------------------
    par_minmax.local_search_flags.validation = false;    
    par_minmax.local_search_flags.inner = false;       % Idem inner (SO) problem
    par_minmax.local_search_flags.outer = false;       % Idem outer (MO) problem

algo_minmax.par_minmax = par_minmax;





% INNER LOOP (maximisation over u)

%--------------------------------------------------------------------------
% algorithm in the form [x,f,exitflag,output] = algo(problem,par_algo)
%--------------------------------------------------------------------------
algo_inner.optimise = @optimise_mpaidea_wrapper;          

%--------------------------------------------------------------------------
% maximum number of function evaluation for the inner loop
%--------------------------------------------------------------------------
    par_mpaidea.nFeValMax = 1000;  
    
%--------------------------------------------------------------------------    
% number of populations, if no adaptive behaviour should set to 1
%-------------------------------------------------------------------------- 
    par_mpaidea.n_populations = 1; 
    
%--------------------------------------------------------------------------    
% number of agents in one population
%--------------------------------------------------------------------------     
    par_mpaidea.n_agents = max(5,problem.dim_u);   
    
% -------------------------------------------------------------------------
% Maximum number of local restart before global restart (for cases when
% only one population is considered and no adaptation of delta_local and
% local/global restart is performed)
% par_mpaidea.max_LR = []; --> adaptation
% -------------------------------------------------------------------------
    par_mpaidea.max_LR = 10; 

% -------------------------------------------------------------------------
% other parameters: default values
% -------------------------------------------------------------------------
    par_mpaidea.plots = 0;                                  % Display plots during run?
    par_mpaidea.text = 0;                                   % Display text during run?
    par_mpaidea.save_pop_DE = 0;                            % Save results of DE to file?
    par_mpaidea.name_save_pop_DE = 'pop_DE_';               % If yes, choose prefix of name for files:
    par_mpaidea.save_local_search = 0;                      % Save results of local search to file?
    par_mpaidea.name_save_local_search = 'minima_fmincon_'; % If yes, choose prefix for name for file:
    par_mpaidea.save_pop_LR = 0;                            % Save populations at local restart (each one saved on a different file)?
    name_save_pop_LR = 'pop_LR_';                           % If yes, choose prefix of name for files:
    par_mpaidea.save_pop_GR = 0;                            % Save populations at global restart (each one saved on a different file)?
    name_save_pop_GR = 'pop_GR_';                           % If yes, choose prefix of name for files:

    par_mpaidea.population=[];                     % initial population, same for all execution. Leave empty to randomize.
    par_mpaidea.DE_strategy = 1;                   % 1: DE/Rand and DE/CurrentToBest; 2: DE/Rand and DE/Best
    par_mpaidea.prob_DE_strategy = 0.5;            % probability or not using DE/Rand
    par_mpaidea.delta_local = 0.1;                 % dimension of the bubble for the local restart of population.
    par_mpaidea.delta_global = 0.1;                % characteristic dimension for the global restart of the population
    par_mpaidea.rho = 0.25;                        % contraction threshold for the population
    par_mpaidea.dd_CRF = 3;                        % parameter for the adaptation of CRF.
    par_mpaidea.F = 0.2;                           % F
    par_mpaidea.CR = 0.8;                          % CR
    par_mpaidea.text = 0;                          % verbosity
    par_mpaidea.record = [1];                      % Fraction(s) of nfevalmax at which the result is saved
    
algo_inner.par = par_mpaidea;



% OUTER LOOP (minimisation over d)

%--------------------------------------------------------------------------
% algorithm in the form [x,f,exitflag,output] = algo(problem,par_algo)
%--------------------------------------------------------------------------
algo_outer.optimise = @optimise_mpaidea_wrapper;    

%--------------------------------------------------------------------------
% maximum number of function evaluation for the outer loop
%--------------------------------------------------------------------------
    par_mpaidea.nFeValMax = 1000;     
        
%--------------------------------------------------------------------------    
% number of populations, if no adaptive behaviour should set to 1
%--------------------------------------------------------------------------       
    par_mpaidea.n_populations = 1;   
        
%--------------------------------------------------------------------------    
% number of agents in one population
%--------------------------------------------------------------------------  
    par_mpaidea.n_agents = max(5,problem.dim_d);   

% -------------------------------------------------------------------------
% Maximum number of local restart before global restart (for cases when
% only one population is considered and no adaptation of delta_local and
% local/global restart is performed)
% par_mpaidea.max_LR = []; --> adaptation
% -------------------------------------------------------------------------
    par_mpaidea.max_LR = 10; 
        
% -------------------------------------------------------------------------
% other parameters: default values
% -------------------------------------------------------------------------    
    par_mpaidea.plots = 0;                                  % Display plots during run?
    par_mpaidea.text = 0;                                   % Display text during run?
    par_mpaidea.save_pop_DE = 0;                            % Save results of DE to file?
    par_mpaidea.name_save_pop_DE = 'pop_DE_';               % If yes, choose prefix of name for files:
    par_mpaidea.save_local_search = 0;                      % Save results of local search to file?
    par_mpaidea.name_save_local_search = 'minima_fmincon_'; % If yes, choose prefix for name for file:
    par_mpaidea.save_pop_LR = 0;                            % Save populations at local restart (each one saved on a different file)?
    name_save_pop_LR = 'pop_LR_';                           % If yes, choose prefix of name for files:
    par_mpaidea.save_pop_GR = 0;                            % Save populations at global restart (each one saved on a different file)?
    name_save_pop_GR = 'pop_GR_';                           % If yes, choose prefix of name for files:

    par_mpaidea.population=[];                % initial population, same for all execution. Leave empty to randomize.
    par_mpaidea.DE_strategy = 1;              % 1: DE/Rand and DE/CurrentToBest; 2: DE/Rand and DE/Best
    par_mpaidea.prob_DE_strategy = 0.5;       % probability or not using DE/Rand
    par_mpaidea.delta_local = 0.1;            % dimension of the bubble for the local restart of population.
    par_mpaidea.delta_global = 0.1;           % characteristic dimension for the global restart of the population
    par_mpaidea.rho = 0.25;                   % contraction threshold for the population
    par_mpaidea.dd_CRF = 3;                   % parameter for the adaptation of CRF.
    par_mpaidea.F = 0.2;                      % F
    par_mpaidea.CR = 0.8;                     % CR
    par_mpaidea.text = 0;                     % verbosity
    par_mpaidea.record = [1];                 % Fraction(s) of nfevalmax at which the result is saved
        
algo_outer.par = par_mpaidea;



% DECOMPOSITION ALGORITHM 

%--------------------------------------------------------------------------
% algorithm in the form [x,f,exitflag,output] = algo(problem,par_algo)
%--------------------------------------------------------------------------
algo_decomposition.optimise = @optimise_mpaidea_wrapper;          

%--------------------------------------------------------------------------
% maximum number of function evaluation for the inner loop
%--------------------------------------------------------------------------
    par_mpaidea.nFeValMax = 1000;  
    
%--------------------------------------------------------------------------    
% number of populations, if no adaptive behaviour should set to 1
%-------------------------------------------------------------------------- 
    par_mpaidea.n_populations = 1; 
    
%--------------------------------------------------------------------------    
% number of agents in one population
%--------------------------------------------------------------------------     
    par_mpaidea.n_agents = max(5,problem.dim_u);   
    
% -------------------------------------------------------------------------
% Maximum number of local restart before global restart (for cases when
% only one population is considered and no adaptation of delta_local and
% local/global restart is performed)
% par_mpaidea.max_LR = []; --> adaptation
% -------------------------------------------------------------------------
    par_mpaidea.max_LR = 10; 

% -------------------------------------------------------------------------
% other parameters: default values
% -------------------------------------------------------------------------
    par_mpaidea.plots = 0;                                  % Display plots during run?
    par_mpaidea.text = 0;                                   % Display text during run?
    par_mpaidea.save_pop_DE = 0;                            % Save results of DE to file?
    par_mpaidea.name_save_pop_DE = 'pop_DE_';               % If yes, choose prefix of name for files:
    par_mpaidea.save_local_search = 0;                      % Save results of local search to file?
    par_mpaidea.name_save_local_search = 'minima_fmincon_'; % If yes, choose prefix for name for file:
    par_mpaidea.save_pop_LR = 0;                            % Save populations at local restart (each one saved on a different file)?
    name_save_pop_LR = 'pop_LR_';                           % If yes, choose prefix of name for files:
    par_mpaidea.save_pop_GR = 0;                            % Save populations at global restart (each one saved on a different file)?
    name_save_pop_GR = 'pop_GR_';                           % If yes, choose prefix of name for files:

    par_mpaidea.population=[];                     % initial population, same for all execution. Leave empty to randomize.
    par_mpaidea.DE_strategy = 1;                   % 1: DE/Rand and DE/CurrentToBest; 2: DE/Rand and DE/Best
    par_mpaidea.prob_DE_strategy = 0.5;            % probability or not using DE/Rand
    par_mpaidea.delta_local = 0.1;                 % dimension of the bubble for the local restart of population.
    par_mpaidea.delta_global = 0.1;                % characteristic dimension for the global restart of the population
    par_mpaidea.rho = 0.25;                        % contraction threshold for the population
    par_mpaidea.dd_CRF = 3;                        % parameter for the adaptation of CRF.
    par_mpaidea.F = 0.2;                           % F
    par_mpaidea.CR = 0.8;                          % CR
    par_mpaidea.text = 0;                          % verbosity
    par_mpaidea.record = [1];                      % Fraction(s) of nfevalmax at which the result is saved
    
algo_decomposition.par = par_mpaidea;



  

%% minmax and/or minmin
[minmin, minmax, LIST, LIST_EXACT] = ebro(problem, algo_minmax, algo_outer, algo_inner, algo_decomposition);


