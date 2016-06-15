clear all; close all; clc
% This is the main file for the Atira problem
%
% Author: Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk

%To do:
%Check agent previous decisions - 2x start?
%Update ramification & agentmovement 1
%Update Atira to change in check & node creation order

%Add the path
addpath(genpath(strcat(pwd,'/AsteroidMainBelt')));
addpath(genpath(strcat(fileparts(fileparts(pwd)),'/Optimisation/AIDMAP')));
addpath(strcat(fileparts(fileparts(pwd)),'/Optimisation'));

%Define the files names for the initialisation
filenames.AsteroidsFileName = 'AsteroidMainBelt/InputFiles/First20kMainBelt.xlsx';
filenames.MatFileName = 'AsteroidMainBelt/InputFiles/MainBeltFirst20kAsteroids.mat';
filenames.NameFile = 'AsteroidMainBelt/InputFiles/MainBeltFirst20kNames.txt';
filenames.epochsnodename = 'AsteroidMainBelt/InputFiles/MainBeltFirst20kEpoch.mat';
filenames.orbitcharsname = 'AsteroidMainBelt/InputFiles/MainBeltFirst20kOrbitChars.mat';

%epoch_start = [10957.5,11016.5,11077.5,11138.5,11200.5,11261.5,];
%epoch_end = [18262.5,18321.5,18382.5,18443.5,18505.5,18566.5];

%epoch_start = [10957.5,10988.5,11016.5,11047.5,11077.5,11108.5,11138.5,11169.5,11200.5,11230.5,11261.5,11291.5];
%epoch_end = [18262.5,18293.5,18321.5,18352.5,18382.5,18413.5,18443.5,18474.5,18505.5,18535.5,18566.5,18596.5];

% Best for 150km
% epoch_start = 11077.5;
% epoch_end = 18382.5;

epoch_start = 10957.5;
epoch_end = 18262.5;

%Define the mean anomalies
startmeananomalies = 0;% :90:270; %270 Best for 150km

%Loop over all the input starting epochs and mean anomalies
for p = 1:length(epoch_start) 
for q = 1:length(startmeananomalies)
    
%Clear the variables obtained during the loop to prevent issues
clearvars -except startmeananomalies p q startorbit filenames epoch_start epoch_end

%Create a diary for this iteration
%diary on
diaryfilename = strcat(['AsteroidMainBelt/Results/First20k/EllipticalStart/MaxdV5/DiaryMainBelt60_M',num2str(startmeananomalies(q)),'Start',strrep(num2str(epoch_start(p)),'.','_'),'_',datestr(now,'yyyymmdd_HHMMSS'),'_','NewRam']);
diary(diaryfilename)

%Initialize the asteroid main belt problem - comment for all asteroid
%mission
InitializeAsteroidsMainBelt(epoch_start(p),epoch_end(p),filenames.AsteroidsFileName,filenames.MatFileName,filenames.NameFile,filenames.epochsnodename,filenames.orbitcharsname);

%Define the starting orbit
%startorbit = [2.85 0	0	0	0	startmeananomalies(q) epoch_start(p)];
startorbit = [1.9512 0.4607	0	0	0	startmeananomalies(q) epoch_start(p)];


%Show the current mean anomaly and start date being evaluated    
disp(char(strcat('Current Mean anomaly:',{' '},num2str(startmeananomalies(q)))));
disp(char(strcat('Current Start Date:',{' '},num2str(epoch_start(p)))));

   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Physarum Options         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options.LinearDilationCoefficient = 5e-3;                       %Linear dilation coefficient 'm'
options.EvaporationCoefficient = 1e-4;                          %Evaporation coefficient 'rho'
options.GrowthFactorVal = 5e-1;                                 %Growth factor 'GF'
options.NumberOfAgents = 5;                                    %Number of virtual agents 'N_agents'
options.RamificationProbability = 0.7;                          %Probability of ramification 'p_ram'
options.RamificationWeight = 1;                                 %Weight on ramification 'lambda'
options.MaximumRadiusRatio = 20;                                %Maximum ratio between the link's radius & the starting radius
options.MinimumRadiusRatio = 1e-3;                              %Maximum ratio between the link's radius & the starting radius
options.StartingRadius = 1;                                     %The starting radius of the veins
options.RamificationAmount = 3;                                 %The number of nodes initially generated for the ramification
options.Generations = 10;                                        %The number of generations
options.Viscosity = 1;                                          %The viscocity of the "fluid" 
options.MinCommonNodesThres = 5;                                %The minimum number of nodes two decision sequences should have in common for a restart to occur
options.IfZeroLength = 1e-15;                                   %Value assigned to the length if it's zero (to prevent flux = inf)
options.MaxChildFindAttempts = 1e6;
options.MinPickProbability = 0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Problem-Specific Options     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options.Targets = textread(filenames.NameFile,'%s')';
options.MaxConsecutiveRes = 0*ones(1, length(options.Targets)); %The maximum number of resonance orbits to each target (set to -1 to ignore)
options.MaxVisits = 1*ones(1, length(options.Targets));           %The maximum nubmer of visists to each target (set to -1 to ignore)                    
options.AttributeIDIndex = [13 12];                             %Index of the attributes that determine the unique ID
options.RootAttrib = [0 startorbit(7)];                                %Attributes of the root  
options.NodeCheckBoundaries = [0.5 0.31 2 365 5];                   %The values used by the MyCreatedNodeCheck file. In this case, it denotes [max dV_dep, min a_per, C for the LT check, max waiting time]  
fitnessfcn = @MyCostFunctionMainBelt;                                   %The function reference to the cost function
options.NodeAttributes = @MyAttributesMainBelt;                         %The class that contains the node attributes
options.MyAttributeCalcFile = @MyAttributeCalcsMainBelt;                %The file that does the additonal calculations wrt the attributes
options.MyNodeIDCheck = @MyNodeCheckMainBelt;                           %The function that checks whether a node can be linked. Can only use the UID
options.MyCreatedNodeCheck = @MyCreatedNodeCheckMainBelt;               %After the node has been found valid using its UID and its structure has been generated, this function checks whether the node itself matches the boundaries
options.MyBestChainFile = @MyBestChainMainBelt;
options.MyEndConditionsFile = @MyEndConditions;
options.EndConditions = {{}};                                           %For use in the MyEndCondtions file
options.RootName = 'Start';

AsteroidsMainBelt = load(filenames.MatFileName);
rootfieldnames = fieldnames(AsteroidsMainBelt.Asteroids.(options.RootName));
for i = 2:length(rootfieldnames)
    AsteroidsMainBelt.Asteroids.(options.RootName).(char(rootfieldnames(i))) = startorbit(i-1);
end

options.AdditonalInputs{1} = AsteroidsMainBelt.Asteroids;
options.AdditonalInputs{2} = 0;                                 %Set to 1 for LT, 0 for HT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Sets input             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tofvalues = 0:40:1.5*365;             %Set the value of the sets. 
sets.tof = mat2cell(ones(length(options.Targets),1)... %Input should be a cell array where each line depicts a target.
   *tofvalues,[ones(length(options.Targets),1)],...    %For this mission, the ToF and the arrival epochs have been used
   [length(tofvalues)]);
load(filenames.epochsnodename)
sets.epochsnode = epochsnode(2:end);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Run the optimiser        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[BestSolution, BestCost, exitflag, output] = optimise_aidmap(fitnessfcn,sets,options);    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Display the result         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PhysarumTreePlot(output.ListNodes)
set(gca,'xcolor','w','ycolor','w','xtick',[],'ytick',[]);

%If no additional nodes have been found, skip to the next starting mean
%anomaly or quit if only one defined
if(length(fieldnames(output.ListNodes))==1)
    %continue
end

%Find solutions with the most asteroids
for i = 1:length(output.Solutions.Nodes);
    asteroidnum(i) = length(output.Solutions.Nodes{i});
end
maxasteroidnumindex = find(asteroidnum==max(asteroidnum));

%Save the solutions with the max number of asteroids found into one array
for i = 1:length(maxasteroidnumindex)
    AllBestSolutions{i,1} = output.Solutions.Nodes{maxasteroidnumindex(i)};
end

%Find the individual costs corresponding to the best solution
for i = 1:length(BestSolution)
    for j = 2:length(BestSolution{i})
      bestnodecosts{i}(j-1) = output.ListNodes.(char(BestSolution{i}(j))).length;
    end
end

%Plot the trajectory
[r] = PlotTrajectories(BestSolution,bestnodecosts,output.ListNodes);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Save the result          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Save only the solution with the most asteroids and the least dV
% filename = strcat(['AsteroidMainBelt/Results/10km/EllipticalStart/MaxdV5/MainBelt60_M',num2str(startmeananomalies(q)),'Startdate',strrep(num2str(epoch_start(p)),'.','_'),num2str(length(AllBestSolutions{1})-1),'Asteroids',num2str(i),'_',num2str(options.NumberOfAgents),'Agents',num2str(options.Generations),'Generations','_',datestr(now,'yyyymmdd_HHMMSS'),'_','NewRam']);
% SaveTrajectorySolution(BestSolution{1},output.ListNodes,strcat(filename));

%Alternative: save all solutions with the max number of asteroids
for i = 1:length(AllBestSolutions)
    filename = strcat(['AsteroidMainBelt/Results/First20k/EllipticalStart/MaxdV5/MainBelt60_M',num2str(startmeananomalies(q)),'Startdate',strrep(num2str(epoch_start(p)),'.','_'),num2str(length(AllBestSolutions{1})-1),'Asteroids',num2str(i),'_',num2str(options.NumberOfAgents),'Agents',num2str(options.Generations),'Generations','_',datestr(now,'yyyymmdd_HHMMSS'),'_','NewRam',num2str(i)]);
    SaveTrajectorySolution(AllBestSolutions{i},output.ListNodes,strcat(filename));
end

%Save the workspace
save(strcat('AsteroidMainBelt/Results/First20k/EllipticalStart/MaxdV5/MainBelt',num2str(options.NumberOfAgents),'Agents',num2str(options.Generations),'Generations','M',num2str(startmeananomalies(q)),'Startdate',strrep(num2str(epoch_start(p)),'.','_'),datestr(now,'yyyymmdd_HHMMSS')));
diary off

end
end