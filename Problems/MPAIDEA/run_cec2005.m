% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%-----------Copyright (C) 2018 University of Strathclyde and Authors-----------
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example of run of optimisation problem of CEC 2005 using MP-AIDEA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
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
% Warm start?
% Provide populations for warm start. The number of populations and number
% of individuals have to be compatible with the number of populations and
% individuals defined above.
% -------------------------------------------------------------------------
options.warm_start = 0;
% Prefix of the names of the file containing the population
name_warm_start = 'pop_GR_';


% -------------------------------------------------------------------------
% Display text during run?
% -------------------------------------------------------------------------
options.text = 1;

% -------------------------------------------------------------------------
% Display plot during run?
% -------------------------------------------------------------------------
options.plots = 1;


% -------------------------------------------------------------------------
% Save results of DE to file?
% All the individuals of all the populations will be saved on a file after
% each generation of the DE
% -------------------------------------------------------------------------
% 1 for yes, 0 for no
options.save_pop_DE = 1;
% If yes, choose the prefix for the name for the file:
options.name_save_pop_DE = 'pop_DE_';


% -------------------------------------------------------------------------
% Save results of local search to file?
% All the local minima are saved to the same file
% -------------------------------------------------------------------------
options.save_local_search = 1;
% If yes, choose the prefix for the name for file:
options.name_save_local_search = 'minima_fmincon_';


% -------------------------------------------------------------------------
% Save populations at local restart (each one saved on a different file)?
% -------------------------------------------------------------------------
options.save_pop_LR = 1;
% If yes, choose prefix of name for files:
name_save_pop_LR = 'pop_LR_';

% -------------------------------------------------------------------------
% Save populations at global restart (each one saved on a different file)?
% -------------------------------------------------------------------------
options.save_pop_GR = 1;
% If yes, choose prefix of name for files:
name_save_pop_GR = 'pop_GR_';

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

% -------------------------------------------------------------------------
% Format for files
% -------------------------------------------------------------------------
% File to save population of DE
if options.save_pop_DE
    options.str = '%8.6e';
    for i = 1 : D
        options.str = [options.str,' ', '%8.6e'];
    end
    options.str = [options.str, '\n'];
    for s = 1 : pop_number
        options.fileID(s) = fopen(strcat(options.name_save_pop_DE, num2str(s),'.txt'),'w');
    end
end

% File to save local searches
if options.save_local_search
    % If yes, uncomment the following and give name to files:
    options.str = '%8.6e';
    for i = 1 : D
        options.str = [options.str,' ', '%8.6e'];
    end
    options.str = [options.str, '\n'];
    for s = 1 : pop_number
        options.fileID2(s) = fopen(strcat(options.name_save_local_search, num2str(s),'.txt'),'w');
    end
end

% File to save local restarts
if options.save_pop_LR
    % If yes, uncomment the following and give name to files:
    options.str2 = '%8.6e';
    for i = 1 : D - 1
        options.str2 = [options.str2,' ', '%8.6e'];
    end
    options.str2 = [options.str2, '\n'];
    for s = 1 : pop_number
        options.fileID3(s) = fopen(strcat(name_save_pop_LR,num2str(s),'.txt'),'w');
    end
end

% File to save global restarts
if options.save_pop_GR
    % If yes, uncomment the following and give name to files:
    options.str2 = '%8.6e';
    for i = 1 : D - 1
        options.str2 = [options.str2,' ', '%8.6e'];
    end
    options.str2 = [options.str2, '\n'];
    for s = 1 : pop_number
        options.fileID4(s) = fopen(strcat(name_save_pop_GR,num2str(s),'.txt'),'w');
    end
end

% -------------------------------------------------------------------------
% Maximum number of function evaluations
% -------------------------------------------------------------------------
options.nFeValMax = nFeValMax;

% Solutions are saved not only when nFeValMax has been reached but also at
% fraction of nFeValMax. Define this fraction in options.record:
options.record = [0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1];

% -------------------------------------------------------------------------
% Populations
% -------------------------------------------------------------------------
% Initialise populations
population = zeros(NP,D,pop_number);

% Population is defined by latin hypercube in all cases, except for warm
% start, when population is defined instead by given files
if options.warm_start == 0
    
    for s = 1 : pop_number
        pop = lhsdesign(NP,D,'criterion','maximin').*repmat(UB-LB,NP,1)+repmat(LB,NP,1);
        population(:,:,s) = pop;
    end
    
else
    
    for s = 1 : pop_number
        population(:,:,s) = textread( strcat( name_warm_start, num2str(1), '.txt') );
    end
    
end

options.population = population;


%% Optimisation

% Function to optimise
fitnessfcn.obj = @(x)benchmark_func(x,func_num);
fitnessfcn.constr = [] ;
fitnessfcn.obj_constr = 0;
fitnessfcn.weighted = 0;

% MP-AIDEA optimisation
[x,fval,exitflag,output] = optimise_mpaidea(fitnessfcn, LB, UB, options);



%% Close file opened for writing 

if options.save_pop_DE
    for s = 1 : pop_number
        fclose(options.fileID(s));
    end
end
if options.save_local_search
    for s = 1 : pop_number
        fclose(options.fileID2(s));
    end
end
if options.save_pop_LR
    for s = 1 : pop_number
        fclose(options.fileID3(s));
    end
end
if options.save_pop_GR
    for s = 1 : pop_number
        fclose(options.fileID4(s));
    end
end

