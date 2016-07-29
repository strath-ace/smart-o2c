% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%-----------Copyright (C) 2016 University of Strathclyde-------------
%
%
%
%% run_atira: This is the main file for the Atira problem. 
% 
% Author(s): Aram Vroom, Marilena Di Carlo and Juan Manuel Romero Martin (2016)
% Email: aram.vroom@strath.ac.uk marilena.di-carlo@strath.ac.uk juan.romero-martin@strath.ac.uk

clear all; close all; clc
rng('shuffle')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Define Paths            % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isunix
    addpath(genpath(strcat(pwd, '/Atira')));
    addpath(genpath(strcat(fileparts(fileparts(pwd)), '/Optimisation/AIDMAP')));
    addpath(strcat(fileparts(fileparts(pwd)), '/Optimisation'));   
    addpath(genpath(strcat(pwd, '/common')));
else
    addpath(genpath(strcat(pwd, '\Atira')));
    addpath(genpath(strcat(fileparts(fileparts(pwd)), '\Optimisation\AIDMAP')));
    addpath(strcat(fileparts(fileparts(pwd)), '\Optimisation'));
    addpath(genpath(strcat(pwd, '\common')));
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Physarum Options         % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options.LinearDilationCoefficient = 1e-3;                       % Linear dilation coefficient 'm' [real number]
options.EvaporationCoefficient = 1e-4;                          % Evaporation coefficient 'rho' [real number]
options.GrowthFactorVal = 5e-3;                                 % Growth factor 'GF' [real number]
options.NumberOfAgents = 20;                                    % Number of virtual agents 'N_agents' [integer]
options.RamificationProbability = 0.7;                          % Probability of ramification 'p_ram' [real number between 0 and 1, where 1 is a 100 probability for an agent to ramificate]
options.RamificationWeight = 1;                                 % Weight on ramification 'lambda' [real number, where a larger value puts more weight on ramification]
options.MaximumRadiusRatio = 2.5;                               % Maximum ratio between the link's radius & the starting radius [real number]
options.MinimumRadiusRatio = 1e-3;                              % Maximum ratio between the link's radius & the starting radius [real number]
options.StartingRadius = 2;                                     % The starting radius of the veins [real number]
options.RamificationAmount = 5;                                 % The number of nodes initially generated for the ramification [integer]
options.Generations = 100;                                      % The number of generations [integer]
options.Viscosity = 1;                                          % The fluid viscocity "mu" [real number]
options.MinCommonNodesThres = 4;                                % The minimum number of nodes two agents in a generation should have in common for a restart to occur [integer]
options.IfZeroLength = 1e-15;                                   % Value assigned to the length if it's zero (to prevent flux = inf) [real number]
options.MaxChildFindAttempts = 1e4;                             % Max number of attempts that will be done to find additional children for a node [integer]      
options.MinPickProbability = 0.1;                               % The minimum probability for a feasible node to be picked before the algorithm changes its method of choosing a child [real number between 0 and 1]
options.GenerateGraphPlot = 0;                                  % Indicator as to whether the algorithm should generate a graph plot animation, where 1 is defined as "yes"
options.GraphPlotFileName = '';                                 % Name of the file that the graph plot animation will be saved as [string]
options.GenerateTreePlot = 0;                                   % Indicator as to whether the algorithm should generate a tree plot, where 1 is defined as "yes"
options.SaveHistory = 0;                                        % Indicator as to whether the algorithm should save the history of the radius of each vein and the path of each agent throughout the simulation, where 1 is defined as "yes"
options.LowMem = 0;                                             % Indicator as to whether the algorithm should use the low-memory version of the searching for new nodes, where 1 is defined as "yes"

SaveDir = 'Atira\IO_Dir\';                                      % Input / Output Directory 


% If the user is using a Linux or Mac version of MATLAB, replace the
% backslashes by forward slashes in the Input / Output directory
if isunix
    SaveDir = strrep(SaveDir,'\','/');                           
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Problem-Specific Options     % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options.Cities = {'neo2003CP20', 'neo2004XZ130', ...            % The list of possible cities [1xC string array, where C is the number of cities]
    'neo1998DK36', 'neo2004JG6', 'neo2005TG45', ...
    'neo2006WE4', 'neo2007EB26', 'neo2008EA32', ...
    'neo2008UL90' , 'neo2010XB11', 'neo2012VE46' , ...
    'neo2013JX28', 'neo2013TQ5', 'neo2014FO47', ...
    'neo2015DR215', 'neo2015ME131'}; 
options.MaxConsecutiveVis = 1*ones(1, length(options.Cities));  % Maximum number of consecutive visits to each city. Set maxima to -1 if no maximum defined [1xC vector of integers, with C being the number of cities]
options.MaxVisits = ones(1, length(options.Cities));            % Maximum number of visits to each city. Set maxima to -1 if no maximum defined [1xC vector of integers, with C being the number of cities]                    
options.AttributeIDIndex = [9 8];                               % Index of the optimisation variables in the MyNodeAttributes class [1xV vector of integer, where V is the number of optimisation variables]
options.RootAttrib = [0 7304.5];                                % Values of the optimisation variables at the root
options.NodeCheckBoundaries = [3 1.5 0.31 2 2*365 5];           % The values used by the MyCreatedNodeCheck file. In this case, it denotes [max dV_dep root, max dV_dep child, min a_per, C for the LT check, max waiting time]
fitnessfcn = @MyCostFunction;                                   % The function reference to the cost function
options.MyNodeAttributes = @MyAttributes;                       % Reference to the file containing the NodeAttributes class, which defines the problems-specific attributes each node has
options.MyAttributeCalcFile = @MyAttributeCalcs;                % The file that does the additonal calculations wrt the attributes
options.MyNodeIDCheck = @MyNodeCheck;                           % The function that checks whether a node can be linked. Can only use the optimisation variables listed in the AttributeIDIndex option
options.MyCreatedNodeCheck = @MyCreatedNodeCheck;               % After the node has been found valid using the optimisation variables and its structure has been generated, this function checks whether the node itself matches the boundaries
options.MyBestChainFile = @MyBestChain;                         % Reference to the file that determines the best chain in a generation
options.EndConditions = {};                                     % End conditions used in the options.MyEndConditionsFile file [cell array]
options.MyEndConditionsFile = @MyEndConditions;                 % Reference to the file that checks whether the end conditions have been reached
options.EndConditions = {{}};                                   % End conditions used in the options.MyEndConditionsFile file [cell array]
options.RootName = 'EARTH';                                     % The name of the root [string]
options.AdditonalInputs = {{}};                                 % Variable that can be used to store any additional information required in one of the user's files [cell array]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Define Sets input        % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tofvalues = 35:10:365;                                
sets.tof = mat2cell(ones(length(options.Cities), 1)... % The structure containing the values that the optimisation
*tofvalues, [ones(length(options.Cities), 1)], ...     % variables can have for each city (where "city" is defined as done 
[length(tofvalues)]);                                  % in the traveling salesman problem). Thus, each field
                                                       % contains a Cx1 cell array, where C is the number of cities and
load('epochsnode.mat')                                 % each cell in turn contains the discrete set of values the optimisation
sets.epochsnode = epochsnode(1:end);                   % variable can have. For this mission, the ToF and the arrival epochs have been used.
													   % To obtain these passing epochs when the asteroid database is changed, AsteroidEpochs.m
													   % should be changed and run first.
   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Run the optimiser        % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run optimise_aidmap
[BestSolution, BestCost, exitflag, output] = optimise_aidmap(fitnessfcn, sets, options);    

% It may be noted that the command window shows the movement of the agents.
% This is shown through the use of the node's unique identifier (UID).
% As only underscores can be used in field names, the UID is made up as follows:
% 3 underscores denote the delimiter between the parent's section of the ID and the child's
% 2 underscores denote the delimiter between the city's name and the optimisation variables
% 1 underscore denotes the delimiter between two optimation variables
% To minimise the length of the UID, the optimisation variables are represented 
% through their index within the sets structure. 
% The first integer after the city's name is the index of the city within the options.Cities array.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Display the result         % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find all solutions with the most asteroids
for i = 1:length(output.Solutions.Nodes);
    asteroidnum(i) = length(output.Solutions.Nodes{i});
end
maxasteroidnumindex = find(asteroidnum==max(asteroidnum));

for i = 1:length(maxasteroidnumindex)
    AllBestSolutions{i, 1} = output.Solutions.Nodes{maxasteroidnumindex(i)};
end

% Plot the solutions with the most asteroids
[r] = PlotTrajectories(AllBestSolutions, output.ListNodes, 1e8, 12, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Save the result          % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save all the solutions that have the max. number of asteroids
for i = 1:length(AllBestSolutions)
    filename = strcat([SaveDir, 'Atira', num2str(length(AllBestSolutions{1})-1), ...
        'Asteroids', num2str(i), '_', num2str(options.NumberOfAgents), 'Agents', ...
        num2str(options.Generations), 'Generations', '_', datestr(now, 'yyyymmdd_HHMMSS')]);
    SaveTrajectorySolution(AllBestSolutions{i}, output.ListNodes, filename);
    ExportSolution(output.ListNodes, AllBestSolutions{i}, filename)
end

% Save the workspace
save(strcat(SaveDir, 'Atira', num2str(options.NumberOfAgents), 'Agents', num2str(options.Generations), 'Generations', '_', datestr(now, 'yyyymmdd_HHMMSS')));

