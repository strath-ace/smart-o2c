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

global initial_flag;
initial_flag = 0;

% Add path to optimiser folder
if isunix
    addpath(genpath('../../Optimisation'))
else
    addpath(genpath('..\..\Optimisation'))
end

% Add path to problem folder
addpath(genpath('CEC2017Constrained'))


%% Choose the CEC 2017 problem 

% Number of the function to optimise. The CEC 2014 competition includes 30
% test functions. func_num must be betwen 1 and 30
func_num_values = 1;

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
% If yes, choose prefix for name for the file:
options.name_save_pop_DE = 'pop_DE_';


% -------------------------------------------------------------------------
% Save results of local search to file?
% All the local minima are saved to the same file
% -------------------------------------------------------------------------
options.save_local_search = 1;
% If yes, choose prefix for name for file:
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



%% Lower and upper boundaries and dimension of the problem


% Lower and upper boundaries of the search space
UB =  100*ones(1,D);
LB = -100*ones(1,D);

% Maximum number of function evaluations
nFeValMax = 20000 * D;


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
        options.fileID3(s) = fopen(strcat(name_save_pop_LR, num2str(s),'.txt'),'w');
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
        options.fileID4(s) = fopen(strcat(name_save_pop_GR, num2str(s),'.txt'),'w');
    end
end

% -------------------------------------------------------------------------
% Maximum number of function evaluations
% -------------------------------------------------------------------------
% Maximum number of function evaluations
options.nFeValMax = nFeValMax;

% Solutions are saved not only when nFeValMax has been reached but also at
% fraction of nFeValMax. Define this fraction in options.record:
options.record = [2000, 1e4, 2e4] * D;

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




%% Example 1: objective and constraints in the same function


% Flag to 1: objective and constraints are in the same function
fitnessfcn.obj_constr = 1;
% How to handle constraints: set to 1 for weighted constraints with fixed
% weights, or to 0 for penalty with no weights
fitnessfcn.weighted = 0;
% If the constraints are handled without weights, then define a tolerance
% for the violation of the equality constraints
fitnessfcn.ceq_eps = 1e-6;
% Weights for penalty
% fitnessfcn.w_ceq = 1000;
% fitnessfcn.w_c = 100;


%% Example 2: objective and constraints are defined in different functions

% % Function to optimise
% fitnessfcn.obj       = @objective;
% % Function of constraints
% fitnessfcn.constr = @constraints;
% % Flag to 0: objective and constraints are NOT in the same function
% fitnessfcn.obj_constr = 0;
% % How to handle constraints: set to 1 for weighted constraints with fixed
% % weights, or to 0 for penalty with no weights
% fitnessfcn.weighted = 0;
% % If the constraints are handled without weights, then define a tolerance
% % for the violation of the equality constraints
% fitnessfcn.ceq_eps = 1e-6;
% % Weights for penalty if fitnessfcn.weighted == 1
% fitnessfcn.w_ceq = 100;
% fitnessfcn.w_c = 100;


%% MP-AIDEA optimisation

for ij = 1 : numel(func_num_values)
    
    % Function to optimise
    fitnessfcn.obj       = @(x)CEC2017(x, func_num_values(ij));
    % Function of constraints
    fitnessfcn.constr = @(x)CEC2017(x, func_num_values(ij));
    
    for n_run = 1 : 25
        
        n_run
        
        [x{n_run},fval{n_run},exitflag,output{n_run}] = optimise_mpaidea(fitnessfcn, LB, UB, options);
        fval{n_run}
        for i = 1 : size(x{n_run},1)
            
            [f(i,n_run),g{n_run}(i,:),h{n_run}(i,:)] = CEC2017(x{n_run}(i,:), func_num_values(ij));
            
        end
        
        % Select only feasible solution
        possible_f = [];
        possible_index = [];
        for i = 1 : 4
            if all(g{n_run}(i,:) <0) && all(abs(h{n_run}(i,:)) < 1e-8)
                possible_index = [possible_index; i];
                possible_f = [possible_f; f(i,n_run)];
            end
        end
        
        % Among the feasible solution, take the one with lower f
        [a,b] = min(possible_f);
        
        % If there are no feasible solutions
        if isempty(possible_f)
            
            % Find the most feasible - how?
            
            % If both equality and inequality constraints are not satisfied
            if all(g{n_run}(i,:) >0) && all(abs(h{n_run}(i,:)) > 1e-8)
                % ?????
                
                % If equality constraints are satisfied but inequality are
                % not satisfied
            elseif all(g{n_run}(i,:) >0) && all(abs(h{n_run}(i,:)) < 1e-8)
                % Choose soltion with lower violation of inequalities - ma
                % sono piu` inequality!
                [~,b] = max(sum(abs(g{n_run}),2));
                a = f(b,n_run);
                A(:,n_run) = [a; g{n_run}(b,:)'; h{n_run}(b,:)'];
                
                % If inequality constraints are satisfied but equality are
                % not satisfied
                
            elseif all(g{n_run}(i,:) <=0) && all(abs(h{n_run}(i,:)) > 1e-8)
                % Choose soltion with lower violation of inequalities - ma
                % sono piu` inequality!
                [~,b] = min(sum(abs(h{n_run}),2));
                a = f(b,n_run);
                A(:,n_run) = [a; g{n_run}(b,:)'; h{n_run}(b,:)'];
            end
            
        else
                A(:,n_run) = [a; g{n_run}(possible_index(b),:)'; h{n_run}(possible_index(b),:)']; 
        end
        
        
        
        
        
    end
    
%     output_file_name = strcat('MP_AIDEA_',num2str(func_num_values(ij)),'_',num2str(D),'_constr.txt');
%     fileID = fopen(output_file_name,'w');
%     % fprintf(fileID, '%12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\r\n', A);
% %     fprintf(fileID, '%12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e \r\n', A');
%     fclose(fileID);
    
    clear f g h A
    
end


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


