clear all; close all; clc
addpath(genpath(strcat(pwd,'/Atira')));
addpath(genpath(strcat(fileparts(fileparts(pwd)),'/Optimisation/AIDMAP')));
addpath(strcat(fileparts(fileparts(pwd)),'/Optimisation'));

% This is the main file for the Atira problem
%
% Author: Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Physarum Options         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options.LinearDilationCoefficient = 5e-3;                       %Linear dilation coefficient 'm'
options.EvaporationCoefficient = 1e-4;                          %Evaporation coefficient 'rho'
options.GrowthFactorVal = 5e-3;                                 %Growth factor 'GF'
options.NumberOfAgents = 2;                                    %Number of virtual agents 'N_agents'
options.RamificationProbability = 0.4;                          %Probability of ramification 'p_ram'
options.RamificationWeight = 1;                                 %Weight on ramification 'lambda'
options.MaximumRadiusRatio = 20;                                %Maximum ratio between the link's radius & the starting radius
options.MinimumRadiusRatio = 1e-3;                              %Maximum ratio between the link's radius & the starting radius
options.StartingRadius = 1;                                     %The starting radius of the veins
options.RamificationAmount = 3;                                 %The number of nodes initially generated for the ramification
options.Generations = 1;                                       %The number of generations
options.Viscosity = 1;                                          %The viscocity of the "fluid" 
options.MinCommonNodesThres = 5;                                %The minimum number of nodes two decision sequences should have in common for a restart to occur
options.IfZeroLength = 1e-15;                                   %Value assigned to the length if it's zero (to prevent flux = inf)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Problem-Specific Options     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options.Targets = {'neo163693', 'neo164294', ...        %Targets the Physarum can choose from
    'neo1998DK36', 'neo2004JG6', 'neo2005TG45',...
    'neo2006WE4', 'neo2007EB26', 'neo2008EA32',...
    'neo2008UL90' ,'neo2010XB11','neo2012VE46' ,...
    'neo2013JX28','neo2013TQ5', 'neo2014FO47', ...
    'neo2015DR215', 'neo2015ME131'}; 
options.MaxConsecutiveRes = 1*ones(1, length(options.Targets)); %The maximum number of resonance orbits to each target (set to -1 to ignore)
options.MaxVisits = ones(1, length(options.Targets));           %The maximum nubmer of visists to each target (set to -1 to ignore)                    
options.AttributeIDIndex = [11 10];                             %Index of the attributes that determine the unique ID
options.RootAttrib = [0 7304.5];                                %Attributes of the root  
options.NodeCheckBoundaries = [3 0.31 2 2*365];                 %The values used by the MyCreatedNodeCheck file.  In this case, it denotes [max dV_dep, min a_per, C for the LT check, max waiting time]
fitnessfcn = @MyCostFunction;                                   %The function reference to the cost function
options.NodeAttributes = @MyAttributes;                         %The class that contains the node attributes
options.MyAttributeCalcFile = @MyAttributeCalcs;                %The file that does the additonal calculations wrt the attributes
options.MyNodeIDCheck = @MyNodeCheck;                           %The function that checks whether a node can be linked. Can only use the UID
options.MyCreatedNodeCheck = @MyCreatedNodeCheck;               %After the node has been found valid using its UID and its structure has been generated, this function checks whether the node itself matches the boundaries
options.MyBestChainFile = @MyBestChain;
options.RootName = 'Earth';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Sets input             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tofvalues = 30:5:365;             %Set the value of the sets. 
sets.tof = mat2cell(ones(16,1)... %Input should be a cell array where each line depicts a target.
   *tofvalues,[ones(16,1)],...    %For this mission, the ToF and the arrival epochs have been used
   [length(tofvalues)]);
load('epochsnode.mat')
sets.epochsnode = epochsnode(2:end);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Run the optimiser        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[BestSolution, BestCost, exitflag, output] = optimise_aidmap(fitnessfcn,sets,options);    
%[output] = optimise_aidmap(fitnessfcn,sets,options);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Display the result         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PhysarumTreePlot(output.ListNodes)
set(gca,'xcolor','w','ycolor','w','xtick',[],'ytick',[]);


%Find solutions with the most asteroids
for i = 1:length(output.Solutions.Nodes);
    asteroidnum(i) = length(output.Solutions.Nodes{i});
end
maxasteroidnumindex = find(asteroidnum==max(asteroidnum));

%Plot the solutions with the most asteroids
for i = 1:length(maxasteroidnumindex)
    AllBestSolutions{i,1} = output.Solutions.Nodes{maxasteroidnumindex(i)};
    %AllBestCosts{i,1} = output.Solutions.Costs{maxasteroidnumindex(i)};
end

%Find the individual costs corresponding to the best solution
for i = 1:length(BestSolution)
    for j = 2:length(BestSolution{i})
      bestnodecosts{i}(j-1) = output.ListNodes.(char(BestSolution{i}(j))).length;
    end
end

[r] = PlotTrajectories(BestSolution,bestnodecosts,output.ListNodes);

%Save all the solutions
for i = 1:length(AllBestSolutions)
    filename = strcat([num2str(length(AllBestSolutions{1})-1),'Asteroids',num2str(i),'_',num2str(options.NumberOfAgents),'Agents',num2str(options.Generations),'Generations','_',datestr(now,'yyyymmdd_HHMMSS'),'_','NewRam']);
    SaveTrajectorySolution(AllBestSolutions{i},output.ListNodes,strcat(filename));
end

%Notes:
%20160510 - GrowthFactor has been changed


