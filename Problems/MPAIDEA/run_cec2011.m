% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%-----------Copyright (C) 2016 University of Strathclyde-------------
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example of run of optimisation problem of CEC 2011 using MP-AIDEA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

% Add path to optimiser folder
if isunix
    addpath(genpath('../../Optimisation'))
else
    addpath(genpath('..\..\Optimisation'))
end

% Add path to problem folder
addpath(genpath('CEC2011'))



%% Choose the CEC 2011 problem 

% Only problem without equality and disequality contraints can be tested using
% MP-AIDEA. The possible problems are: 1, 2, 3, 5, 6, 7, 10, 12, 13

% Number of the function to optimise. Uncomment one of the following lines:
func_num = 1;
% func_num = 2;
% func_num = 3;
% func_num = 5;
% func_num = 6;
% func_num = 7;
% func_num = 10;
% func_num = 12;
% func_num = 13;


%% Set the parameters for MP-AIDEA

% -------------------------------------------------------------------------
% Number of populations
% -------------------------------------------------------------------------
pop_number = 4;

% -------------------------------------------------------------------------
% Number of individuals in each population
% -------------------------------------------------------------------------
NP = 30;

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
% Maximum number of local restart before global restart (for cases when 1
% population is considered and no adaptation of delta_local and
% local/global restart is performed)
% -------------------------------------------------------------------------
options.max_LR = [];
% options.max_LR = 5;

% -------------------------------------------------------------------------
% Choose the DE strategies. 
% -------------------------------------------------------------------------
% The DE evolution of MP-AIDEA uses two DE strategies, with probability
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
% Example: if options_DE_strategy was set to 1, options.prob_DE_strategy
% defines the probability of using DE/Rand rather than DE/CurrentToBest
options.prob_DE_strategy = 0.5;


% -------------------------------------------------------------------------
% Parameter for the adaptation of CR and F
% -------------------------------------------------------------------------
% Value of CR (if empty, MP-AIDEA-ALR adapt will adapt it during the process)
options.CR = 0.5;
% options.CR = [];

% Value of F (if empty, MP-AIDEA will adapt it during the process)
options.F = 0.5;
% options.F = [];

% If options.CR and options.F are empty, define CRF for adaptation of CR
% and F:
options.dd_CRF = 3;


% -------------------------------------------------------------------------
% Display text during run?
% -------------------------------------------------------------------------
options.text = 1;


%% Lower and upper boundaries and dimension of the problem

% Dimension of the problem - The dimension of the problem and the lower and
% upper boundaries depend on the given problem
switch func_num
    
    case 1
        D = 6;
        LB = -6.4 * ones(1,D);
        UB = 6.35 * ones(1,D);
        
    case 2
        D = 30;
        UB = [4 4 pi ];
        LB = [0 0 0];
        for coord = 4 : D
            UB(coord) = 4 + 0.25 * (coord - 4)/3;
            LB(coord) = - 4 - 0.25 * (coord - 4)/3;
        end
        
    case 3
        D = 1;
        LB = 0.6 * ones(1,D);
        UB = 0.9 * ones(1,D);
        
    case 5
        D = 30;
        UB = [4 4 pi ];
        LB = [0 0 0];
        
        for coord = 4 : D
            UB(coord) = 4 + 0.25 * (coord - 4)/3;
            LB(coord) = - 4 - 0.25 * (coord - 4)/3;
        end
        
    case 6
        D = 30;
        
        UB = [4 4 pi ];
        LB = [0 0 0];
        
        for coord = 4 : D
            UB(coord) = 4 + 0.25 * (coord - 4)/3;
            LB(coord) = - 4 - 0.25 * (coord - 4)/3;
        end
        
    case 7
        D = 20;
        
        LB = 0 * ones(1,D);
        UB = 2*pi * ones(1,D);
        
    case 10
        D = 12;
        
        LB = [0.2*ones(1,D/2) -180*ones(1,D/2)];
        UB = [1*ones(1,D/2) 180*ones(1,D/2)];
                
    case 12
        D = 26;
        
        LB = [1900 2.5  0 0 100 100 100 100 100 100 0.01 0.01 0.01 0.01 0.01 0.01 ...
              1.1 1.1 1.05 1.05 1.05 -pi -pi -pi -pi -pi];
        UB = [2300 4.05 1 1 500 500 500 500 500 600 0.99 0.99 0.99 0.99 0.99 0.99 ...
              6    6   6     6     6  pi  pi  pi  pi  pi];

        
    case 13
        D = 22;
        
        LB = [-1000 3 0 0 100 100 30 400 800 0.01 0.01 0.01 0.01 0.01 1.05 1.05 1.15 1.7 -pi -pi -pi -pi];
        UB = [0 5 1 1 400 500 300 1600 2200 0.9 0.9 0.9 0.9 0.9 6 6 6.5 291 pi pi pi pi];

end


% Maximum number of function evaluations
nFeValMax = 150000;


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
switch func_num
    case {1, 2, 3, 5, 6, 7}
        fitnessfcn = @(x)bench_func(x,func_num);
    case 10
        fitnessfcn = @(x)antennafunccircular(x,[50,120],180,0.5);
    case 12 
        load messengerfull.mat
        fitnessfcn = @(x)messengerfull(x, MGADSMproblem);
    case 13
        load cassini2.mat
        fitnessfcn = @(x)cassini2(x,MGADSMproblem);
end

% MP-AIDEA optimisation
[x,fval,exitflag,output] = optimise_mpaidea(fitnessfcn, LB, UB, options);




