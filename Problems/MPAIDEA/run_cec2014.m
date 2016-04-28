%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example of run of optimisation problem of CEC 2014 using MP-AIDEA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

% Add path to optimiser folder
addpath(genpath('..\..\Optimisation'))
addpath(genpath('CEC2014'))

% Mexify cpp file
mex CEC2014/cec14_func.cpp 

%% Choose the CEC 2014 problem

% Number of the function to optimise. The CEC 2014 competition includes 30
% test functions. func_num must be betwen 1 and 30
func_num = 1;

% Dimension of the problem - Choose between 10, 30, 50 and 100 dimensions
D = 10;


%% Set the parameters for MP-AIDEA

% -------------------------------------------------------------------------
% Number of populations
% -------------------------------------------------------------------------
pop_number = 4;

% -------------------------------------------------------------------------
% Number of individuals in each population
% -------------------------------------------------------------------------
NP = 10;

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
% Choose the DE strategies. 
% -------------------------------------------------------------------------
% The DE evolution of MP-AIDEA uses two DE strategies, with probability
% defined by options.prob_DE_strategy.
% DE/Rand, DE/CurrentToBest and DE/Best are well know DE strategies.
% DE/Archive is ..... (explain)
% Uncomment the following line for DE strategies DE/Rand and DE/CurrentToBest:
options.DE_strategy = 1;
% Uncomment the following line for DE strategies DE/Rand and DE/Best:
% options.DE_strategy = 2;
% Uncomment the following line for DE strategies DE/Best and DE/Archive
% options.DE_strategy = 3;
% Uncomment the following line for DE strategies DE/Rand and DE/Archive
% options.DE_strategy = 4;


% -------------------------------------------------------------------------
% Probability of having DE strategy 1 rather than DE strategy 2 (DE
% strategies 1 and 2 have been defined in options.DE_strategy)
% -------------------------------------------------------------------------
% Example: if options_DE_strategy was set to 1, options.prob_DE_strategy
% defines the probability of using DE/Rand rather than DE/CurrentToBest
options.prob_DE_strategy = 0.5;


% -------------------------------------------------------------------------
% Parameter for the adaptation of CRF
% -------------------------------------------------------------------------
options.dd_CRF = 3;

% 
options.plot_flag = 0;

% 
options.text = 0;

%% CEC 2014 guidelines

% Lower and upper boundaries of the search space
UB =  100*ones(1,D);
LB = -100*ones(1,D);

% Maximum number of function evaluations
nFeValMax = 10000 * D;


%% MP-AIDEA inputs

% Maximum number of function evaluations
options.nFeValMax = nFeValMax;





% Initialise populations
population = zeros(NP,D,pop_number);

for s = 1 : pop_number
    pop = lhsdesign(NP,D,'criterion','maximin').*repmat(UB-LB,NP,1)+repmat(LB,NP,1);
    population(:,:,s) = pop;
end


options.population = population;


%% Optimisation

fitnessfcn = @(x)cec14_func(x,func_num);

[x,fval,exitflag,output] = optimise_mpaidea(fitnessfcn, D, LB, UB, options);




