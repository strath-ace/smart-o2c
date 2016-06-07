clear all; close all; clc
% This is the main file for the Atira problem
%
% Author: Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk

%To do:
%Check agent previous decisions - 2x start?

%Add the path
addpath(genpath(strcat(pwd,'/AsteroidMainBelt')));
addpath(genpath(strcat(fileparts(fileparts(pwd)),'/Optimisation/AIDMAP')));
addpath(strcat(fileparts(fileparts(pwd)),'/Optimisation'));

%Define the files names for the initialisation
filenames.AsteroidsFileName = 'AsteroidMainBelt/DiameterGreater150.xlsx';
filenames.MatFileName = 'AsteroidMainBelt/MainBelt60Asteroids.mat';
filenames.NameFile = 'AsteroidMainBelt/MainBelt60Names.txt';
filenames.epochsnodename = 'AsteroidMainBelt/MainBelt60Epoch.mat';
filenames.orbitcharsname = 'AsteroidMainBelt/MainBelt60OrbitChars.mat';

%Define starting & end epochs
epoch_start = [10957.5, 11016.5, 11077.5, 11138.5,11200.5, 11261.5];
epoch_end = [18262.5,18321.5, 18382.5, 18443.5,18505.5, 18566.5];

%Define the mean anomalies
startmeananomalies = 0:45:270;

%Loop over all the input starting epochs and mean anomalies
for p = 1:length(epoch_start) 
for q = 1:length(startmeananomalies)
    
%Clear the variables obtained during the loop to prevent issues
clearvars -except startmeananomalies p q startorbit filenames epoch_start epoch_end

%Create a diary for this iteration
diary on
diaryfilename = strcat(['AsteroidMainBelt/Results/Multi-Start-and-M-30y-set/DiaryMainBelt60_M',num2str(startmeananomalies(q)),'Start',strrep(num2str(epoch_start(p)),'.','_'),'_',datestr(now,'yyyymmdd_HHMMSS'),'_','NewRam']);
diary(diaryfilename)

%Initialize the asteroid main belt problem
InitializeAsteroidsMainBelt(epoch_start(p),epoch_end(p),filenames.AsteroidsFileName,filenames.MatFileName,filenames.NameFile,filenames.epochsnodename,filenames.orbitcharsname);

%Define the starting orbit
startorbit = [2.85	0	0	0	0	startmeananomalies(q) epoch_start(p)];

%Show the current mean anomaly and start date being evaluated    
disp(char(strcat('Current Mean anomaly:',{' '},num2str(startmeananomalies(q)))));
disp(char(strcat('Current Start Date:',{' '},num2str(epoch_start(p)))));

   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Physarum Options         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options.LinearDilationCoefficient = 5e-3;                       %Linear dilation coefficient 'm'
options.EvaporationCoefficient = 1e-4;                          %Evaporation coefficient 'rho'
options.GrowthFactorVal = 5e-1;                                 %Growth factor 'GF'
options.NumberOfAgents = ;                                    %Number of virtual agents 'N_agents'
options.RamificationProbability = 0.7;                          %Probability of ramification 'p_ram'
options.RamificationWeight = 1;                                 %Weight on ramification 'lambda'
options.MaximumRadiusRatio = 20;                                %Maximum ratio between the link's radius & the starting radius
options.MinimumRadiusRatio = 1e-3;                              %Maximum ratio between the link's radius & the starting radius
options.StartingRadius = 1;                                     %The starting radius of the veins
options.RamificationAmount = 3;                                 %The number of nodes initially generated for the ramification
options.Generations = 1;                                       %The number of generations
options.Viscosity = 1;                                          %The viscocity of the "fluid" 
options.MinCommonNodesThres = 7;                                %The minimum number of nodes two decision sequences should have in common for a restart to occur
options.IfZeroLength = 1e-15;                                   %Value assigned to the length if it's zero (to prevent flux = inf)
options.MaxChildFindAttempts = 1e5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Problem-Specific Options     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options.Targets = textread(filenames.NameFile,'%s')';
options.MaxConsecutiveRes = 0*ones(1, length(options.Targets)); %The maximum number of resonance orbits to each target (set to -1 to ignore)
options.MaxVisits = 1*ones(1, length(options.Targets));           %The maximum nubmer of visists to each target (set to -1 to ignore)                    
options.AttributeIDIndex = [12 11];                             %Index of the attributes that determine the unique ID
options.RootAttrib = [0 startorbit(7)];                                %Attributes of the root  
options.NodeCheckBoundaries = [3 0.31 2 365];                   %The values used by the MyCreatedNodeCheck file. In this case, it denotes [max dV_dep, min a_per, C for the LT check, max waiting time]  
fitnessfcn = @MyCostFunctionMainBelt;                                   %The function reference to the cost function
options.NodeAttributes = @MyAttributesMainBelt;                         %The class that contains the node attributes
options.MyAttributeCalcFile = @MyAttributeCalcsMainBelt;                %The file that does the additonal calculations wrt the attributes
options.MyNodeIDCheck = @MyNodeCheckMainBelt;                           %The function that checks whether a node can be linked. Can only use the UID
options.MyCreatedNodeCheck = @MyCreatedNodeCheckMainBelt;               %After the node has been found valid using its UID and its structure has been generated, this function checks whether the node itself matches the boundaries
options.MyBestChainFile = @MyBestChainMainBelt;
options.EndTarget = {};
options.RootName = 'Start';

AsteroidsMainBelt = load(filenames.MatFileName);
rootfieldnames = fieldnames(AsteroidsMainBelt.Asteroids.(options.RootName));
for i = 2:length(rootfieldnames)
    AsteroidsMainBelt.Asteroids.(options.RootName).(char(rootfieldnames(i))) = startorbit(i-1);
end

options.AdditonalInputs{1} = AsteroidsMainBelt.Asteroids;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Sets input             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tofvalues = 200:50:610;             %Set the value of the sets. 
sets.tof = mat2cell(ones(length(options.Targets),1)... %Input should be a cell array where each line depicts a target.
   *tofvalues,[ones(length(options.Targets),1)],...    %For this mission, the ToF and the arrival epochs have been used
   [length(tofvalues)]);
load(filenames.epochsnodename)
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

if(length(fieldnames(output.ListNodes))==1)
    continue
end

%Find solutions with the most asteroids
for i = 1:length(output.Solutions.Nodes);
    asteroidnum(i) = length(output.Solutions.Nodes{i});
end
maxasteroidnumindex = find(asteroidnum==max(asteroidnum));

%Plot the solutions with the most asteroids
for i = 1:length(maxasteroidnumindex)
    AllBestSolutions{i,1} = output.Solutions.Nodes{maxasteroidnumindex(i)};
end

%Find the individual costs corresponding to the best solution
for i = 1:length(BestSolution)
    for j = 2:length(BestSolution{i})
      bestnodecosts{i}(j-1) = output.ListNodes.(char(BestSolution{i}(j))).length;
    end
end

[r] = PlotTrajectories(BestSolution,bestnodecosts,output.ListNodes);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Save the result          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for i = 1:length(BestSolution)
    filename = strcat(['AsteroidMainBelt/Results/Multi-Start-and-M-30y-set/MainBelt60_M',num2str(startmeananomalies(q)),'Startdate',strrep(num2str(epoch_start(p)),'.','_'),num2str(length(AllBestSolutions{1})-1),'Asteroids',num2str(i),'_',num2str(options.NumberOfAgents),'Agents',num2str(options.Generations),'Generations','_',datestr(now,'yyyymmdd_HHMMSS'),'_','NewRam']);
    SaveTrajectorySolution(BestSolution{1},output.ListNodes,strcat(filename));
    %SaveTrajectorySolution(AllBestSolutions{i},output.ListNodes,strcat(filename));
%end

save(strcat('AsteroidMainBelt/Results/Multi-Start-and-M-30y-set/MainBelt',num2str(options.NumberOfAgents),'Agents',num2str(options.Generations),'Generations','M',num2str(startmeananomalies(q)),'Startdate',strrep(num2str(epoch_start(p)),'.','_'),datestr(now,'yyyymmdd_HHMMSS')));
%diary off

end
end