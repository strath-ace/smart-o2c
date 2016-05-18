clear all; close all; clc
addpath(genpath(strcat(pwd,'/AsteroidMainBelt')));
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
options.NumberOfAgents = 1;                                    %Number of virtual agents 'N_agents'
options.RamificationProbability = 0.7;                          %Probability of ramification 'p_ram'
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
options.Targets = {'Ceres','Pallas','Juno','Vesta','Astraea','Hebe','Iris','Flora','Metis','Hygiea','Parthenope','Victoria','Egeria','Irene','Eunomia','Psyche','Melpomene','Fortuna','Massalia','Kalliope','Thalia','Themis','Bellona','Amphitrite','Urania','Euphrosyne','Circe','Leukothea','Atalante','Fides','Leda','Laetitia','Harmonia','Daphne','Isis','Eugenia','Hestia','Aglaja','Doris','Pales','Nemausa','Europa','Kalypso','Alexandra','Melete','Mnemosyne','Elpis','Ausonia','Leto','Hesperia','Panopaea','Galatea','Diana','Terpsichore','Io','Semele','Thisbe','Julia','Antiope','Aegina','Undina','Minerva','Aurora','Arethusa','Aegle','Ianthe','Klymene','Artemis','Dione','Ate','Lomia','Lachesis','Johanna','Nemesis','Antigone','Elektra','Sophrosyne','Meliboea','Juewa','Siwa','Lumen','Vibilia','Adeona','Lucina','Protogeneia','Nuwa','Bertha','Xanthippe','Aemilia','Eva','Loreley','Ophelia','Ino','Andromache','Iduna','Eucharis','Eunike','Lamberta','Kolga','Nausikaa','Prokne','Philomela','Dynamene','Pompeja','Hersilia','Dido','Isolda','Medea','Kleopatra','Eos','Athamantis','Asterope','Hypatia','Vanadis','Germania','Eukrate','Aletheia','Aline','Adorea','Sapientia','Adelheid','Emma','Polyxo','Bamberga','Gudrun','Desiderata','Hermentaria','Dembowska','Ornamenta','Eleonora','Liguria','Ninina','Carlova','Corduba','Palma','Ursula','Siegena','Aquitania','Charybdis','Thia','Aspasia','Chloris','Aurelia','Diotima','Hippo','Gyptis','Patientia','Argentina','Papagena','Hedwig','Kreusa','Comacina','Veritas','Cava','Marion','Princetonia','Davida','Armida','Brixia','Herculina','Scheila','Marianna','Elfriede','Notburga','Zelinda','Gerlinde','Wratislavia','Alauda','Interamnia','Erminia','Boliviana','Mandeville','Winchester','Faina','Pulcova','Tatjana','Tanete','Berbericia','Hohensteina','Ani','Hispania','Tauris','Freda'}; 
options.MaxConsecutiveRes = 0*ones(1, length(options.Targets)); %The maximum number of resonance orbits to each target (set to -1 to ignore)
options.MaxVisits = ones(1, length(options.Targets));           %The maximum nubmer of visists to each target (set to -1 to ignore)                    
options.AttributeIDIndex = [11 10];                             %Index of the attributes that determine the unique ID
options.RootAttrib = [0 10957.5];                                %Attributes of the root  
options.NodeCheckBoundaries = [3 2 2 2*365];               %The values used by the MyCreatedNodeCheck file. In this case, it denotes [max dV_dep, min a_per, C for the LT check, max waiting time]  
fitnessfcn = @MyCostFunctionMainBelt;                                   %The function reference to the cost function
options.NodeAttributes = @MyAttributesMainBelt;                         %The class that contains the node attributes
options.MyAttributeCalcFile = @MyAttributeCalcsMainBelt;                %The file that does the additonal calculations wrt the attributes
options.MyNodeIDCheck = @MyNodeCheckMainBelt;                           %The function that checks whether a node can be linked. Can only use the UID
options.MyCreatedNodeCheck = @MyCreatedNodeCheckMainBelt;               %After the node has been found valid using its UID and its structure has been generated, this function checks whether the node itself matches the boundaries
options.MyBestChainFile = @MyBestChainMainBelt;
options.EndTarget = {};
options.RootName = 'Start';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Sets input             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tofvalues = 180:50:2*365;             %Set the value of the sets. 
sets.tof = mat2cell(ones(length(options.Targets),1)... %Input should be a cell array where each line depicts a target.
   *tofvalues,[ones(length(options.Targets),1)],...    %For this mission, the ToF and the arrival epochs have been used
   [length(tofvalues)]);
load('AsteroidMainBelt/epochsnode.mat')
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Save the result          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(AllBestSolutions)
    filename = strcat(['MainBelt_',num2str(length(AllBestSolutions{1})-1),'Asteroids',num2str(i),'_',num2str(options.NumberOfAgents),'Agents',num2str(options.Generations),'Generations','_',datestr(now,'yyyymmdd_HHMMSS'),'_','NewRam']);
    SaveTrajectorySolution(AllBestSolutions{i},output.ListNodes,strcat(filename));
end

%Notes:
%20160510 - GrowthFactor has been changed


