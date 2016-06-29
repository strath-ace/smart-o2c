function [InitializedInputs, ListNodes] = InitializePhysarum(fitnessfcn,options,sets)
% This is the main file of the physarum solver. 
% It contains the logic for the algorithm and the solver parameters.
%
% Inputs:
% * fitnessfcn : Reference to the fitness function
% * options    : The structure containg the options set by the user
% * sets       : The structure with the possible values for the attributes
%                put in the UID
%
% Outputs: 
% * InitializedInputs     : Structure containing the PhysarumSolver inputs
% * ListNodes             : The structure containing all the nodes
%
% Author: Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk

disp('Initializing Physarum...')

%Solver Parameters:
InitializedInputs = struct(...
                'LinearDilationCoefficient',          options.LinearDilationCoefficient,  ...   %Linear dilation coefficient 'm'
                'EvaporationCoefficient',             options.EvaporationCoefficient,  ...      %Evaporation coefficient 'rho'
                'GrowthFactor',                       options.GrowthFactorVal,  ...             %Growth factor 'GF'
                'NumberOfAgents',                     options.NumberOfAgents,  ...              %Number of virtual agents 'N_agents'
                'RamificationProbability',            options.RamificationProbability, ...      %Probability of ramification 'p_ram'
                'RamificationWeight',                 options.RamificationWeight,  ...          %Weight on ramification 'lambda'
                'MaximumRadiusRatio',                 options.MaximumRadiusRatio,  ...          %Maximum radius of the veins
                'MinimumRadiusRatio',                 options.MinimumRadiusRatio,  ...          %Minimum radius of the veins
                'StartingRadius',                     options.StartingRadius,  ...              %The starting radius of the veins
                'RamificationAmount',                 options.RamificationAmount,  ...          %The number of nodes initially generated for the ramification
                'PossibleDecisions',                  {options.Targets},  ...                   %The list of possible targets
                'MaxConsecutiveRes',                  options.MaxConsecutiveRes,  ...           %Maximum number of consecutive resonance orbits to each target
                'MaxVisits',                          options.MaxVisits,  ...                   %Maximum number of visits to each target
                'RootAttrib',                         options.RootAttrib,   ...                 %Attributes of the root
                'NodeCheckBoundaries',                options.NodeCheckBoundaries  , ...        %The values used by the MyCreatedNodeCheck file
                'Generations',                        options.Generations, ...                  %The number of generations
                'Viscosity',                          options.Viscosity, ...                    %The viscocity of the "fluid" 
                'MinCommonNodesThres',                options.MinCommonNodesThres,  ...         %The minimum number of nodes two decision sequences should have in common for a restart to occur
                'IfZeroLength',                       options.IfZeroLength, ...                 %Value assigned to the length if it's zero (to prevent flux = inf)
                'MaxChildFindAttempts',               options.MaxChildFindAttempts, ...         %Max number of attempts that will be done to find additional children for a node
                'CostFunction',                       fitnessfcn,  ...                          %Reference to the fitness function file
                'NodeAttributes',                     options.NodeAttributes, ...               %Reference to the file containing the NodeAttributes class
                'MyAttributeCalcFile',                options.MyAttributeCalcFile, ...          %Reference to the file performing the attribute calculations
                'NodeIDCheckFile',                    options.MyNodeIDCheck, ...                %Reference to the file that checks the feasibility of a node's child using solely the unique ID
                'CreatedNodeCheckFile',               options.MyCreatedNodeCheck, ...           %Reference to the file that checks the feasibility of a node's child after the child's structure has been created
                'BestChainFile',                      options.MyBestChainFile, ...              %Reference to the file that determines the best chain
                'EndConditionsFile',                  options.MyEndConditionsFile, ...          %Reference toe the file that checks whether the end conditions have been reached
                'EndConditions',                      {options.EndConditions}, ...              %End conditions used in the options.MyEndConditionsFile file
                'AttributeIDIndex',                   options.AttributeIDIndex, ...             %Index of the attributes that determine the unique ID
                'AdditionalInputs',                   {options.AdditonalInputs}, ...            %Made into cell in case multiple additional outpu
                'Sets',                               {sets}, ...                               %The sets of the possible attributes incorporated in the ID
                'RootName',                           options.RootName, ...                     %The name of the root
                'MinPickProbability',                 options.MinPickProbability, ...            %The minimum probability for a feasible node to be picked before the algorithm changes its method of choosing a child
                'GenerateGraphPlot',                  options.GenerateGraphPlot, ...
                'SaveHistory',                        options.SaveHistory ...
            );
        
      
%%%Error Checking%%%

%Display error if the vectors with number of possible decisions, max. number of consecutive
%resonance orbits & the max. number of visists is not equal
if ~(length(InitializedInputs.MaxConsecutiveRes)==length(InitializedInputs.PossibleDecisions))
    error('Check size of PossibleDecisions, MaxConsecutiveRes and MaxVisits')
end

%Check if the number of boundaries and stepsizes correspond
if ~(length(InitializedInputs.AttributeIDIndex)==length(fieldnames(sets)))
    error('Check size of AttributeIDIndex and sets')
end

%Add nodes that can be selected to the Inputs structure
InitializedInputs.PossibleListNodes = PossChilds(options.Targets,sets);

%Create the list of nodes
ListNodes = CreateListNodes(InitializedInputs);

clearvars -except ListNodes InitializedInputs
disp('Initialization Complete.')

end

