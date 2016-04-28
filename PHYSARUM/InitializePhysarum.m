function [InitializedInputs, ListNodes] = InitializePhysarum()
% This is the main file of the physarum solver. 
% It contains the logic for the algorithm and the solver parameters.
%
% Inputs:
% * 
%
% Outputs: 
% * Inputs         : Structure containing the PhysarumSolver inputs
%
% Author: Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk

disp('Initializing Physarum...')

%Retrieve the user settings
PhysarumSettings

%Solver Parameters:
InitializedInputs = struct('LowThrust',               LowThrust,  ... %Set to 1 for low-thrust, 0 for high-thrust
                'LinearDilationCoefficient',          LinearDilationCoefficient,  ... %Linear dilation coefficient 'm'
                'EvaporationCoefficient',             EvaporationCoefficient,  ... %Evaporation coefficient 'rho'
                'GrowthFactor',                       GrowthFactorVal,  ... %Growth factor 'GF'
                'NumberOfAgents',                     NumberOfAgents,  ... %Number of virtual agents 'N_agents'
                'RamificationProbability',            RamificationProbability, ... %Probability of ramification 'p_ram'
                'RamificationWeight',                 RamificationWeight,  ... %Weight on ramification 'lambda'
                'MaximumRadiusRatio',                 MaximumRadiusRatio,  ... %Maximum radius of the veins
                'MinimumRadiusRatio',                 MinimumRadiusRatio,  ... %Minimum radius of the veins
                'StartingRadius',                     StartingRadius,  ... %The starting radius of the veins
                'RamificationAmount',                 RamificationAmount,  ... %The number of nodes initially generated for the ramification
                'PossibleDecisions',                  {Targets},  ... %The list of possible targets
                'MaxConsecutiveRes',                  {MaxConsecutiveRes},  ... %Maximum number of consecutive resonance orbits to each target
                'MaxVisits',                          {MaxVisits},  ... %Maximum number of visits to each target
                'RootChar',                           RootChar,   ... %Characteristic of the root
                'Generations',                        Generations, ... %The number of generations
                'Viscosity',                          Viscosity, ... %The viscocity of the "fluid" 
                'DeterminingCharacteristic',          DeterminingCharacteristic, ... %The index of the determining characteristic in the 'characteristics' field
                'MinCommonNodesThres',                MinCommonNodesThres,  ... %The minimum number of nodes two decision sequences should have in common for a restart to occur
                'IfZeroLength',                       IfZeroLength, ... %Value assigned to the length if it's zero (to prevent flux = inf)
                'CostFunction',                       CostFunction,  ...
                'NodeAttributes',                     NodeAttributes, ...
                'ProjectDirectory',                   ProjectDirectory, ...
                'AttributeIDIndex',                   AttributeIDIndex ... %Index of the attributes that determine the unique ID
        );      
        
%Include Project Directory
addpath(InitializedInputs.ProjectDirectory);

%%%Error Checking%%%

%Display error if the vectors with number of possible decisions, max. number of consecutive
%resonance orbits & the max. number of visists is not equal
if ~(length(InitializedInputs.MaxConsecutiveRes)==length(InitializedInputs.PossibleDecisions))
    error('Check size of PossibleDecisions, MaxConsecutiveRes and MaxVisits')
end

%Check if the number of boundaries and stepsizes correspond
if ~(length(InitializedInputs.AttributeIDIndex)==length(charvalues))
    error('Check size of AttributeIDIndex and charvalues')
end

%Add nodes that can be selected to the Inputs structure
InitializedInputs.PossibleListNodes = PossNodes(Targets,charvalues);

%Create the list of nodes
ListNodes = CreateListNodes(InitializedInputs);

clearvars -except ListNodes InitializedInputs
disp('Initialization Complete.')

end

