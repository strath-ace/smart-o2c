% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%-----------Copyright (C) 2016 University of Strathclyde-------------
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example of run of optimisation problem of CEC 2005 using MP-AIDEA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc
global nFeVal
% Add path to optimiser folder
if isunix
    addpath(genpath('../../Optimisation'))
else
    addpath(genpath('..\..\Optimisation'))
end

% Add path to problem folder
addpath(genpath('CEC2005'))

% Add global variable required by CEC 2005 test functions
global initial_flag
initial_flag = 0;

%% Choose the CEC 2005 problem

% Number of the function to optimise. The CEC 2005 competition includes 25
% test functions. func_num must be betwen 1 and 25
func_num = 1;

% Dimension of the problem - Choose between 10, 30 and 50 dimensions
D = 10;


%% Set the parameters for MP-AIDEA

% -------------------------------------------------------------------------
% Number of populations
% -------------------------------------------------------------------------
pop_number = 4;

% -------------------------------------------------------------------------
% Number of individuals in each population
% -------------------------------------------------------------------------
NP = D;

% -------------------------------------------------------------------------
% Dimension of the bubble for the local restart. If empty, MP-AIDEA will
% adapt it
% -------------------------------------------------------------------------
options.delta_local = [];
% options.delta_local = 0.1;

% -------------------------------------------------------------------------
% Distance from cluster of global minima for global restart of the
% population
% -------------------------------------------------------------------------
options.delta_global = 0.1;

% -------------------------------------------------------------------------
% Threshold for contraction of the population
% -------------------------------------------------------------------------
options.rho = 0.2;

% -------------------------------------------------------------------------
% Maximum number of local restart before global restart (for cases when
% only one population is considered and no adaptation of delta_local and
% local/global restart is performed)
% -------------------------------------------------------------------------
options.max_LR =[];
% options.max_LR = 5;

% -------------------------------------------------------------------------
% Choose the Differential Evolution (DE) strategies. 
% -------------------------------------------------------------------------
% The DE of MP-AIDEA uses two DE strategies, with probability
% defined by options.prob_DE_strategy (see later).
% DE/Rand, DE/CurrentToBest and DE/Best are well know DE strategies.
% Uncomment the following line for DE strategies DE/Rand and DE/CurrentToBest:
options.DE_strategy = 1;
% Uncomment the following line for DE strategies DE/Rand and DE/Best:
% options.DE_strategy = 2;


% -------------------------------------------------------------------------
% Probability of having DE strategy 1 rather than DE strategy 2 (DE
% strategies 1 and 2 have been defined in options.DE_strategy)
% -------------------------------------------------------------------------
% Example: if options.DE_strategy was set to 1, options.prob_DE_strategy
% defines the probability of using DE/Rand rather than DE/CurrentToBest
options.prob_DE_strategy = 0.5;


% -------------------------------------------------------------------------
% Parameter for the adaptation of CR and F
% -------------------------------------------------------------------------
% Value of CR (if empty, MP-AIDEA-ALR adapt will adapt it during the process)
% options.CR = 0.5;
options.CR = [];

% Value of F (if empty, MP-AIDEA will adapt it during the process)
% options.F = 0.5;
options.F = [];

% If options.CR and options.F are empty, define CRF for adaptation of CR
% and F:
options.dd_CRF = 3;


% -------------------------------------------------------------------------
% Display text during run?
% -------------------------------------------------------------------------
options.text = 0;



%% Lower and upper boundaries of the search space

% The lower and upper boundaries of the search space depend on the given problem
switch func_num
    
    case {1, 2, 3, 4, 5, 6, 14}
        UB =  100*ones(1,D);
        LB = -100*ones(1,D);
        
    case 7
        UB =  600*ones(1,D);
        LB =  zeros(1,D);
        options.no_bounds = 1;
        
    case 8 
        UB =  32*ones(1,D);
        LB = -32*ones(1,D);
        
    case {9, 10, 13, 15, 16, 17 18, 19, 20, 21, 22, 23, 24}
        UB =  5*ones(1,D);
        LB = -5*ones(1,D);
        
    case 11
        UB =  0.5*ones(1,D);
        LB = -0.5*ones(1,D);
        
    case 12
        UB =  pi*ones(1,D);
        LB = -pi*ones(1,D);
        
    case 25
        UB =  2*ones(1,D);
        LB =  5*ones(1,D);      
end


% Maximum number of function evaluations
nFeValMax = 10000 * D;


%% MP-AIDEA inputs

% Maximum number of function evaluations
options.nFeValMax = nFeValMax;

% Solutions are saved not only when nFeValMax has been reached but also at
% fraction of nFeValMax. Define this fraction in options.record:
options.record = [0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1];

% Initialise populations
population = zeros(NP,D,pop_number);

for s = 1 : pop_number
    pop = lhsdesign(NP,D,'criterion','maximin').*repmat(UB-LB,NP,1)+repmat(LB,NP,1);
    population(:,:,s) = pop;
end


options.population = population;


%% Optimisation

% Function to optimise
fitnessfcn.obj = @(x)benchmark_func(x,func_num);
fitnessfcn.constr = [] ;
fitnessfcn.obj_constr = 0;

% MP-AIDEA optimisation
[x,fval,exitflag,output] = optimise_mpaidea(fitnessfcn, LB, UB, options);




