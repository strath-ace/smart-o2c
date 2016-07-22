%% run_atira: This is the main file for the Atira problem
%
%% Author: Aram Vroom (2016)
% Email:  aram.vroom@strath.ac.uk

clear all; close all; clc
rng('shuffle')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Define Paths            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isunix
    addpath(genpath(strcat(pwd,'/Atira')));
    addpath(genpath(strcat(fileparts(fileparts(pwd)),'/Optimisation/AIDMAP')));
    addpath(strcat(fileparts(fileparts(pwd)),'/Optimisation'));    
else
    addpath(genpath(strcat(pwd,'\Atira')));
    addpath(genpath(strcat(fileparts(fileparts(pwd)),'\Optimisation\AIDMAP')));
    addpath(strcat(fileparts(fileparts(pwd)),'\Optimisation'));
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Physarum Options         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options.LinearDilationCoefficient = 1e-3;                       %Linear dilation coefficient 'm' [real number]
options.EvaporationCoefficient = 1e-4;                          %Evaporation coefficient 'rho' [real number]
options.GrowthFactorVal = 5e-3;                                 %Growth factor 'GF' [real number]
options.NumberOfAgents = 2;                                     %Number of virtual agents 'N_agents' [integer]
options.RamificationProbability = 0.7;                          %Probability of ramification 'p_ram' [real number between 0 and 1, where 1 is a 100 probability for an agent to ramificate]
options.RamificationWeight = 1;                                 %Weight on ramification 'lambda' [real number, where a larger value puts more weight on ramification]
options.MaximumRadiusRatio = 2.5;                               %Maximum ratio between the link's radius & the starting radius [real number]
options.MinimumRadiusRatio = 1e-3;                              %Maximum ratio between the link's radius & the starting radius [real number]
options.StartingRadius = 2;                                     %The starting radius of the veins [real number]
options.RamificationAmount = 5;                                 %The number of nodes initially generated for the ramification [integer]
options.Generations = 1;                                        %The number of generations [integer]
options.Viscosity = 1;                                          %The fluid viscocity "mu" [real number]
options.MinCommonNodesThres = 4;                                %The minimum number of nodes two agents in a generation should have in common for a restart to occur [integer]
options.IfZeroLength = 1e-15;                                   %Value assigned to the length if it's zero (to prevent flux = inf) [real number]
options.MaxChildFindAttempts = 1e4;                             %Max number of attempts that will be done to find additional children for a node [integer]      
options.MinPickProbability = 0.1;                               %The minimum probability for a feasible node to be picked before the algorithm changes its method of choosing a child [real number between 0 and 1]
options.GenerateGraphPlot = 0;                                  %Indicator as to whether the algorithm should generate a graph plot animation, where 1 is defined as "yes"
options.GraphPlotFileName = '';                                 %Name of the file that the graph plot animation will be saved as [string]
options.GenerateTreePlot = 0;                                   %Indicator as to whether the algorithm should generate a tree plot, where 1 is defined as "yes"
options.SaveHistory = 0;                                        %Indicator as to whether the algorithm should save the history of the radius of each vein and the path of each agent throughout the simulation, where 1 is defined as "yes"

if isunix
    SaveDir = 'Atira/IO_Dir/';                                  %Output Directory    
else
    SaveDir = 'Atira\IO_Dir\';                                  %Output Directory  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Problem-Specific Options     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options.Cities = {'neo2003CP20', 'neo2004XZ130', ...            %The list of possible cities [1xC string array, where C is the number of cities]
    'neo1998DK36', 'neo2004JG6', 'neo2005TG45',...
    'neo2006WE4', 'neo2007EB26', 'neo2008EA32',...
    'neo2008UL90' ,'neo2010XB11','neo2012VE46' ,...
    'neo2013JX28','neo2013TQ5', 'neo2014FO47', ...
    'neo2015DR215', 'neo2015ME131'}; 
options.MaxConsecutiveVis = 1*ones(1, length(options.Cities));  %Maximum number of consecutive visits to each city. Set maxima to -1 if no maximum defined [1xC vector of integers, with C being the number of cities]
options.MaxVisits = ones(1, length(options.Cities));            %Maximum number of visits to each city. Set maxima to -1 if no maximum defined [1xC vector of integers, with C being the number of cities]                    
options.AttributeIDIndex = [9 8];                               %Index of the optimisation variables in the MyNodeAttributes class [1xV vector of integer, where V is the number of optimisation variables]
options.RootAttrib = [0 7304.5];                                %Values of the optimisation variables at the root
options.NodeCheckBoundaries = [3 1.5 0.31 2 2*365 5];           %The values used by the MyCreatedNodeCheck file. In this case, it denotes [max dV_dep root, max dV_dep child, min a_per, C for the LT check, max waiting time]
fitnessfcn = @MyCostFunction;                                   %The function reference to the cost function
options.MyNodeAttributes = @MyAttributes;                       %Reference to the file containing the NodeAttributes class, which defines the problems-specific attributes each node has
options.MyAttributeCalcFile = @MyAttributeCalcs;                %The file that does the additonal calculations wrt the attributes
options.MyNodeIDCheck = @MyNodeCheck;                           %The function that checks whether a node can be linked. Can only use the UID
options.MyCreatedNodeCheck = @MyCreatedNodeCheck;               %After the node has been found valid using its UID and its structure has been generated, this function checks whether the node itself matches the boundaries
options.MyBestChainFile = @MyBestChain;                         %Reference to the file that determines the best chain
options.EndConditions = {};                                     %End conditions used in the options.MyEndConditionsFile file [cell array]
options.MyEndConditionsFile = @MyEndConditions;                 %Reference to the file that checks whether the end conditions have been reached
options.EndConditions = {{}};                                   %End conditions used in the options.MyEndConditionsFile file [cell array]
options.RootName = 'EARTH';                                     %The name of the root [string]
options.AdditonalInputs = {{}};                                 %Variable that can be used to store any additional information required in one of the user's files [cell array]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Define Sets input        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tofvalues = 35:10:365;                                
sets.tof = mat2cell(ones(length(options.Cities),1)... %The structure containing the values that the optimisation
*tofvalues,[ones(length(options.Cities),1)],...       %variables can have for each city (where "city" is defined as done 
[length(tofvalues)]);                                 %in the traveling salesman problem). Thus, each field
                                                      %contains a Cx1 cell array, where C is the number of cities and
load('epochsnode.mat')                                %each cell in turn contains the discrete set of values the optimisation
sets.epochsnode = epochsnode(1:end);                  %variable can have
                                                      %For this mission, the ToF and the arrival epochs have been used
   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Run the optimiser        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Run optimise_aidmap
[BestSolution, BestCost, exitflag, output] = optimise_aidmap(fitnessfcn,sets,options);    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Display the result         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Find all solutions with the most asteroids
for i = 1:length(output.Solutions.Nodes);
    asteroidnum(i) = length(output.Solutions.Nodes{i});
end
maxasteroidnumindex = find(asteroidnum==max(asteroidnum));

for i = 1:length(maxasteroidnumindex)
    AllBestSolutions{i,1} = output.Solutions.Nodes{maxasteroidnumindex(i)};
end

%Plot the solutions with the most asteroids
[r] = PlotTrajectories(AllBestSolutions,output.ListNodes);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Save the result          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Save all the solutions that have the max. number of asteroids
for i = 1:length(AllBestSolutions)
    filename = strcat([SaveDir,'Atira',num2str(length(AllBestSolutions{1})-1),...
        'Asteroids',num2str(i),'_',num2str(options.NumberOfAgents),'Agents',...
        num2str(options.Generations),'Generations','_',datestr(now,'yyyymmdd_HHMMSS')]);
    SaveTrajectorySolution(AllBestSolutions{i},output.ListNodes,filename);
    ExportSolution(output.ListNodes, AllBestSolutions{i}, filename)
end

%Save the workspace
save(strcat(SaveDir,'Atira',num2str(options.NumberOfAgents),'Agents',num2str(options.Generations),'Generations','_',datestr(now,'yyyymmdd_HHMMSS')));




