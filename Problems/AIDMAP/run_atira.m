clear all; close all; clc
addpath(genpath(fileparts(fileparts(fileparts(pwd)))));
% This is the main file for the Atira problem
%
% Author: Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Sets input             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tofvalues = 30:10:365;             %Set the value of the sets. 
sets.tof = mat2cell(ones(12,1)...  %Input should be a cell array where each line depicts a target.
    *tofvalues,[ones(12,1)],...
    [length(tofvalues)]);
load('epochsnode.mat')
sets.epochsnode = epochsnode;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Options             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options.Targets = {'neo163693', 'neo164294', ...        %Targets the Physarum can choose from
    'neo1998DK36', 'neo2004JG6', 'neo2005TG45',...
    'neo2006WE4', 'neo2007EB26', 'neo2008EA32',...
    'neo2008UL90' ,'neo2010XB11','neo2012VE46' ,...
    'neo2013JX28'}; 
options.MaxConsecutiveRes = 0*ones(1, length(options.Targets));%The maximum number of resonance orbits to each target (set to -1 to ignore)
options.MaxVisits = ones(1, length(options.Targets));        %The maximum nubmer of visists to each target (set to -1 to ignore)                    
options.AttributeIDIndex = [11 10];                     %Index of the attributes that determine the unique ID
options.LowThrust = 1;                                  %Low Thrust Flag. Set to 1 for low-thrust, 0 for high-thrust
options.LinearDilationCoefficient = 20;                 %Linear dilation coefficient 'm'
options.EvaporationCoefficient = 0;                     %Evaporation coefficient 'rho'
options.GrowthFactorVal = 0;                            %Growth factor 'GF'
options.NumberOfAgents = 3;                             %Number of virtual agents 'N_agents'
options.RamificationProbability = 0.4;                 %Probability of ramification 'p_ram'
options.RamificationWeight = 1;                         %Weight on ramification 'lambda'
options.MaximumRadiusRatio = 1000;                      %Maximum ratio between the link's radius & the starting radius
options.MinimumRadiusRatio = 1e-3;                      %Maximum ratio between the link's radius & the starting radius
options.StartingRadius = 1;                             %The starting radius of the veins
options.RamificationAmount = 5;                         %The number of nodes initially generated for the ramification
options.RootAttrib = [0 0];                             %Attributes of the root  
options.Generations = 1;                                %The number of generations
options.Viscosity = 1;                                  %The viscocity of the "fluid" 
options.DeterminingAttribute = 1;                       %The index of the determining attribute in the 'attributes' field
options.MinCommonNodesThres = 7;                        %The minimum number of nodes two decision sequences should have in common for a restart to occur
options.IfZeroLength = 1e-15;                           %Value assigned to the length if it's zero (to prevent flux = inf)
fitnessfcn = @MyCostFunction;                           %The function reference to the cost function
options.NodeAttributes = @MyAttributes;                 %The class that contains the node attributes
options.MyAttributeCalcFile = @MyAttributeCalcs;     %The file that does the additonal calculations wrt the attributes
options.ProjectDirectory = 'C:\Users\ckb16114\Desktop\Internship\Code\Developing\smart-o2c\Problems\AIDMAP'; %The project directory

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Run the optimiser        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[output] = optimise_aidmap(fitnessfcn,sets,options);
%[x,fval,exitflag,output] = optimise_aidmap(fitnessfcn,sets,options)    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Display the result         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PhysarumTreePlot(output.ListNodes)
set(gca,'xcolor','w','ycolor','w','xtick',[],'ytick',[]);




%To do: use Parent kep_trans to find r_dep & v_dep
