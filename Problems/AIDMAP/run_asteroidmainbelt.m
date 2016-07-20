clearvars -except k; close all; clc
rng('shuffle')
% This is the main file for the Atira problem
%
% Author: Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk

%To do:
%Update ramification & agentmovement 1
%Redo ALL solutions (due to bugs)
%Fix if multiple best solutions found in 1 generation

%Add the path and the files that contain the data on the asteroids
if isunix
    addpath(genpath(strcat(pwd,'/AsteroidMainBelt')));
    addpath(genpath(strcat(fileparts(fileparts(pwd)),'/Optimisation/AIDMAP')));
    addpath(strcat(fileparts(fileparts(pwd)),'/Optimisation'));
    
    filenames.AsteroidsFileName = 'AsteroidMainBelt/InputFiles/37Asteroids.xlsx';
    filenames.MatFileName = 'AsteroidMainBelt/InputFiles/37Asteroids.mat';
    filenames.NameFile = 'AsteroidMainBelt/InputFiles/37AsteroidsNames.txt';
    filenames.epochsnodename = 'AsteroidMainBelt/InputFiles/37AsteroidsEpochs.mat';
    filenames.orbitcharsname = 'AsteroidMainBelt/InputFiles/37AsteroidsOrbitChars.mat';
else
    addpath(genpath(strcat(pwd,'\AsteroidMainBelt')));
    addpath(genpath(strcat(fileparts(fileparts(pwd)),'\Optimisation\AIDMAP')));
    addpath(strcat(fileparts(fileparts(pwd)),'\Optimisation'));
    
    filenames.AsteroidsFileName = 'AsteroidMainBelt\InputFiles\37Asteroids.xlsx';
    filenames.MatFileName = 'AsteroidMainBelt\InputFiles\37Asteroids.mat';
    filenames.NameFile = 'AsteroidMainBelt\InputFiles\37AsteroidsNames.txt';
    filenames.epochsnodename = 'AsteroidMainBelt\InputFiles\37AsteroidsEpochs.mat';
    filenames.orbitcharsname = 'AsteroidMainBelt\InputFiles\37AsteroidsOrbitChars.mat';
end

if isunix
    SaveDir = 'AsteroidMainBelt/Results/10km/NewStartDate/EdelBaum/'; %NOT USED
else
    SaveDir = 'AsteroidMainBelt\Results\10km\NewStartDate\EdelBaum\'; %NOT USED
end

epoch_start = 10594.35;
epoch_end = 17894.35;

%Define the starting orbit
startorbit = [2.0533 0.5130	0	0	150	339.35 epoch_start];

%Initialize the asteroid main belt problem: retrieves asteroid nodal
%passing epochs and saves them to a file
InitializeAsteroidsMainBelt(epoch_start,epoch_end,filenames.AsteroidsFileName,filenames.MatFileName,filenames.NameFile,filenames.epochsnodename,filenames.orbitcharsname);
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Physarum Options         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options.LinearDilationCoefficient = 5e-3;                       %Linear dilation coefficient 'm'
options.EvaporationCoefficient = 1e-3;                          %Evaporation coefficient 'rho'
options.GrowthFactorVal = 5e-3;                                 %Growth factor 'GF'
options.NumberOfAgents = 1;                                    %Number of virtual agents 'N_agents'
options.RamificationProbability = 0.7;                          %Probability of ramification 'p_ram'
options.RamificationWeight = 1;                                 %Weight on ramification 'lambda'
options.MaximumRadiusRatio = 2.5;                                %Maximum ratio between the link's radius & the starting radius
options.MinimumRadiusRatio = 1e-3;                              %Maximum ratio between the link's radius & the starting radius
options.StartingRadius = 2;                                     %The starting radius of the veins
options.RamificationAmount = 5;                                 %The number of nodes initially generated for the ramification
options.Generations = 1;                                        %The number of generations
options.Viscosity = 1;                                          %The viscocity of the "fluid" 
options.MinCommonNodesThres = 5;                                %The minimum number of nodes two decision sequences should have in common for a restart to occur
options.IfZeroLength = 1e-15;                                   %Value assigned to the length if it's zero (to prevent flux = inf)
options.MaxChildFindAttempts = 2e4; %2.5e4;
options.MinPickProbability = 0.05;
options.GenerateGraphPlot = 0;
options.SaveHistory = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Problem-Specific Options     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options.Targets = textread(filenames.NameFile,'%s')';
options.MaxConsecutiveRes = 0*ones(1, length(options.Targets)); %The maximum number of resonance orbits to each target (set to -1 to ignore)
options.MaxVisits = 1*ones(1, length(options.Targets));           %The maximum nubmer of visists to each target (set to -1 to ignore)                    
options.AttributeIDIndex = [13 12];                             %Index of the attributes that determine the unique ID
options.RootAttrib = [0 startorbit(7)];                                %Attributes of the root  
options.NodeCheckBoundaries = [0.75 0.31 2 2*365 5];                   %The values used by the MyCreatedNodeCheck file. In this case, it denotes [max dV_dep, min a_per, C for the LT check, max waiting time]  
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Sets input             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tofvalues = 4:10:2*365;                                %Set the value of the sets. 
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

if (exitflag~=0)
    if (options.GenerateGraphPlot ~= 0)
        PhysarumGraphPlot(options, output.ListNodes, output.History,'Test5');
    end
        
    %Plot the trajectory
    [r] = PlotTrajectories(BestSolution,output.ListNodes);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Save the result          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Save the solution with the most asteroids and the least dV
    filename = strcat([SaveDir,'MainBelt_M',strrep(num2str(startorbit(6)),'.','_'),'Startdate',strrep(num2str(epoch_start),'.','_'),num2str(length(BestSolution{1})-1),'Asteroids',num2str(i),'_',num2str(options.NumberOfAgents),'Agents',num2str(options.Generations),'Generations','_',datestr(now,'yyyymmdd_HHMMSS'),'_','NewRam']);
    SaveTrajectorySolution(BestSolution{1},output.ListNodes,filename);

end

%Remove unnecessary outputs
nodenames = fieldnames(output.ListNodes);
for i = 1:length(nodenames)
    output.ListNodes.(char(nodenames(i))) = rmfield(output.ListNodes.(char(nodenames(i))),'ChildValidityTracker');
end  

%Save the workspace
save(strcat(SaveDir,'MainBelt',num2str(options.NumberOfAgents),'Agents',num2str(options.Generations),'Generations','M',strrep(num2str(startorbit(6)),'.','_'),'Startdate',strrep(num2str(epoch_start),'.','_'),datestr(now,'yyyymmdd_HHMMSS')));
