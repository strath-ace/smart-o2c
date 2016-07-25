% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%-----------Copyright (C) 2016 University of Strathclyde-------------
%
% This is the main file for the Main Belt problem
% 
% Author(s): Aram Vroom, Marilena Di Carlo and Juan Manuel Romero Martin (2016)
% Email: aram.vroom@strath.ac.uk marilena.di-carlo@strath.ac.uk juan.romero-martin@strath.ac.uk
%
%% References:
% * JPL Small-Body Database Search Engine
% http://ssd.jpl.nasa.gov/sbdb_query.cgi#x

clear all; close all; clc
rng('shuffle')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Define Paths           % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add the paths
if isunix
    addpath(genpath(strcat(pwd, '/AsteroidMainBelt')));
    addpath(genpath(strcat(fileparts(fileparts(pwd)), '/Optimisation/AIDMAP')));
    addpath(strcat(fileparts(fileparts(pwd)), '/Optimisation'));
    addpath(genpath(strcat(pwd, '/common')));
else
    addpath(genpath(strcat(pwd, '\AsteroidMainBelt')));
    addpath(genpath(strcat(fileparts(fileparts(pwd)), '\Optimisation\AIDMAP')));
    addpath(strcat(fileparts(fileparts(pwd)), '\Optimisation'));    
    addpath(genpath(strcat(pwd, '\common')));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Physarum Options         % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options.LinearDilationCoefficient = 5e-3;                       % Linear dilation coefficient 'm' [real number]
options.EvaporationCoefficient = 1e-4;                          % Evaporation coefficient 'rho' [real number]
options.GrowthFactorVal = 5e-3;                                 % Growth factor 'GF' [real number]
options.NumberOfAgents = 20;                                    % Number of virtual agents 'N_agents' [integer]
options.RamificationProbability = 0.7;                          % Probability of ramification 'p_ram' [real number between 0 and 1, where 1 is a 100 probability for an agent to ramificate]
options.RamificationWeight = 1;                                 % Weight on ramification 'lambda' [real number, where a larger value puts more weight on ramification]
options.MaximumRadiusRatio = 2.5;                               % Maximum ratio between the link's radius & the starting radius [real number]
options.MinimumRadiusRatio = 1e-3;                              % Maximum ratio between the link's radius & the starting radius [real number]
options.StartingRadius = 2;                                     % The starting radius of the veins [real number]
options.RamificationAmount = 5;                                 % The number of nodes initially generated for the ramification [integer]
options.Generations = 40;                                       % The number of generations [integer]
options.Viscosity = 1;                                          % The fluid viscocity "mu" [real number]
options.MinCommonNodesThres = 5;                                % The minimum number of nodes two agents in a generation should have in common for a restart to occur [integer]
options.IfZeroLength = 1e-15;                                   % Value assigned to the length if it's zero (to prevent flux = inf) [real number]
options.MaxChildFindAttempts = 2e4;                             % Max number of attempts that will be done to find additional children for a node [integer]      
options.MinPickProbability = 0.05;                              % The minimum probability for a feasible node to be picked before the algorithm changes its method of choosing a child [real number between 0 and 1]
options.GenerateGraphPlot = 0;                                  % Indicator as to whether the algorithm should generate a graph plot animation, where 1 is defined as "yes"
options.GraphPlotFileName = '';                                 % Name of the file that the graph plot animation will be saved as [string]
options.GenerateTreePlot = 0;                                   % Indicator as to whether the algorithm should generate a tree plot, where 1 is defined as "yes"
options.SaveHistory = 0;                                        % Indicator as to whether the algorithm should save the history of the radius of each vein and the path of each agent throughout the simulation, where 1 is defined as "yes"
options.LowMem = 0;                                             % Indicator as to whether the algorithm should use the low-memory version of searching for new nodes, where 1 is defined as "yes". Using the low-memory version is slower, but requires less memory. 

SaveDir = 'AsteroidMainBelt\IO_Dir\';                           % Input / Output Directory   


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Problem-Specific Options     % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The name of the XLS file that contains the orbital elements of the
% asteroids. Each row should be built up as follows:
% [Name a e i OM W M0 t0],  where a = Semimajor axis [AU], 
% e = Eccentricity, i = Inclination [deg], OM = Asc. Node/raan [deg], 
% W = Arg. Perigee [deg], M0 = Mean anomoly, M at time given t0 [deg] 
% t0 = Time at which M0 is given [MJD]
filenames.AsteroidsFileName = 'AsteroidMainBelt\IO_Dir\DiameterGreater10km_Reduced.xlsx';

% The name of the output files that will hold the data on the asteroids as
% calculated by the InitialiseMainBelt script
filenames.MatFileName = 'AsteroidMainBelt\IO_Dir\DiameterGreater10km_Reduced.mat';
filenames.NameFile = 'AsteroidMainBelt\IO_Dir\DiameterGreater10km_Reduced_Names.txt';
filenames.epochsnodename = 'AsteroidMainBelt\IO_Dir\DiameterGreater10km_Reduced_Epochs.mat';


% If the user is using a Linux or Mac version of MATLAB, replace the
% backslashes by forward slashes in the filenames and Input / Output
% directory
if isunix    
    filenames.AsteroidsFileName = strrep(filenames.AsteroidsFileName,'\','/');
    filenames.MatFileName = strrep(filenames.MatFileName,'\','/');
    filenames.NameFile = strrep(filenames.NameFile,'\','/');
    filenames.epochsnodename = strrep(filenames.epochsnodename,'\','/');    
    SaveDir = strrep(SaveDir,'\','/');     
end

% Start & end epoch of the mission [MJD2000]
epoch_start = 10594.35;
epoch_end = 17894.35;

% Define the starting orbit [a (AU) e (-) i (deg) OM (deg) W (deg) M0 (deg)
% t0 (MJD2000)]. Note that t0 here is in MJD2000, while the XLS uses MJD
startorbit = [2.0533 0.5130	0	0	150	339.35 epoch_start];

% Initialise the asteroid main belt problem: retrieves asteroid nodal
% passing epochs and asteroid names, and saves them to a file
InitialiseAsteroidsMainBelt(epoch_start, epoch_end, filenames);
  
options.Cities = textread(filenames.NameFile, '%s')';               % The list of possible cities [1xC string array, where C is the number of cities]
options.MaxConsecutiveVis = 0*ones(1, length(options.Cities));      % Maximum number of consecutive visits to each city. Set maxima to -1 if no maximum defined [1xC vector of integers, with C being the number of cities]
options.MaxVisits = 1*ones(1, length(options.Cities));              % Maximum number of visits to each city. Set maxima to -1 if no maximum defined [1xC vector of integers, with C being the number of cities]                    
options.AttributeIDIndex = [11 10];                                 % Index of the optimisation variables in the MyNodeAttributes class [1xV vector of integer, where V is the number of optimisation variables]
options.RootAttrib = [0 startorbit(7)];                             % Values of the optimisation variables at the root
options.NodeCheckBoundaries = [0.75 0.31 2 2*365 5];                % The values used by the MyCreatedNodeCheckMainBelt file. In this case, it denotes [max dV_dep, min a_per, C for the LT check, max waiting time]  
fitnessfcn = @MyCostFunctionMainBelt;                               % The function reference to the cost function
options.MyNodeAttributes = @MyAttributesMainBelt;                   % Reference to the file containing the NodeAttributes class, which defines the problems-specific attributes each node has
options.MyAttributeCalcFile = @MyAttributeCalcsMainBelt;            % The file that does the additonal calculations wrt the attributes
options.MyNodeIDCheck = @MyNodeCheckMainBelt;                       % The function that checks whether a node can be linked. Can only use the optimisation variables listed in the AttributeIDIndex option
options.MyCreatedNodeCheck = @MyCreatedNodeCheckMainBelt;           % After the node has been found valid using the optimisation variables and its structure has been generated, this function checks whether the node itself matches the boundaries
options.MyBestChainFile = @MyBestChainMainBelt;                     % Reference to the file that determines the best chain in a generation
options.MyEndConditionsFile = @MyEndConditionsMainBelt;             % Reference to the file that checks whether the end conditions have been reached
options.EndConditions = {{}};                                       % End conditions used in the options.MyEndConditionsFile file [cell array]
options.RootName = 'Start';                                         % The name of the root [string]


% Load the asteroid structure and add the starting location to it
AsteroidsMainBelt = load(filenames.MatFileName);
AsteroidsMainBelt.Asteroids.(options.RootName) = CelestialBody('Start', startorbit(1), startorbit(2), startorbit(3), startorbit(4), startorbit(5), startorbit(6), startorbit(7));


% Save the structure in the AdditionalInputs cell array, such that it can be
% used in other files 
options.AdditonalInputs{1} = AsteroidsMainBelt.Asteroids;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Sets input             % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tofvalues = 4:10:2*365;                                 % The structure containing the values that the optimisation
sets.tof = mat2cell(ones(length(options.Cities), 1)...  % variables can have for each city (where "city" is defined as done 
   *tofvalues, [ones(length(options.Cities), 1)], ...   % in the traveling salesman problem). Thus, each field
   [length(tofvalues)]);                                % contains a Cx1 cell array, where C is the number of cities and
load(filenames.epochsnodename)                          % each cell in turn contains the discrete set of values the optimisation
sets.epochsnode = epochsnode;                           % variable can have. For this mission, the ToF and the arrival epochs have been used


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Run the optimiser        % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run optimise_aidmap
[BestSolution, BestCost, exitflag, output] = optimise_aidmap(fitnessfcn, sets, options);    

% It may be noted that the command window shows the movement of the agents.
% This is shown through the use of the node's unique identifier (UID).
% As only underscores can be used in field names, the UID is made up as follows:
% 3 underscores seperate the parent's section of the ID and the child's
% 2 underscores define the difference between the city's name and the optimisation variables
% 1 underscores define the difference between two optimation variables
% To minimise the length of the UID, the optimisation variables are represented 
% through their index within the sets structure. 
% The first integer after the city's name is the index of the city within the options.Cities array.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Display the result         % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
% Plot the trajectory
[r] = PlotTrajectories(BestSolution, output.ListNodes, 1e8, 12, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Save the result          % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save the solution with the most asteroids and the least dV
filename = strcat([SaveDir, 'MainBelt_', 'Startdate', strrep(num2str(epoch_start), '.', '_'), '_', num2str(length(BestSolution)-1), 'Asteroids', '_', num2str(options.NumberOfAgents), 'Agents', num2str(options.Generations), 'Generations', '_', datestr(now, 'yyyymmdd_HHMMSS')]);
SaveTrajectorySolution(BestSolution, output.ListNodes, filename);
ExportSolution(output.ListNodes, BestSolution, filename)

% Remove unnecessary outputs to reduce workspace .mat file size
nodenames = fieldnames(output.ListNodes);
for i = 1:length(nodenames)
    output.ListNodes.(char(nodenames(i))) = rmfield(output.ListNodes.(char(nodenames(i))), 'ChildValidityTracker');
end  

% Save the workspace
save(strcat(SaveDir, 'MainBelt', num2str(options.NumberOfAgents), 'Agents', num2str(options.Generations), 'Generations', 'M', strrep(num2str(startorbit(6)), '.', '_'), 'Startdate', strrep(num2str(epoch_start), '.', '_'), '_', datestr(now, 'yyyymmdd_HHMMSS')));
