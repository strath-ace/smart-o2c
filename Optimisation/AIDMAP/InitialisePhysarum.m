function [InitialisedInputs, ListNodes] = InitialisePhysarum(fitnessfcn, options, sets)
%% InitialisePhysarum: This algorithm prepares the workspace needed for the AIDMAP solver
% 
%% Inputs:
% * fitnessfcn : Reference to the fitness function
% * options    : The structure containg the options set by the user
% * sets       : The structure with the possible values for the attributes
%                put in the UID
% 
%% Outputs: 
% * InitialisedInputs     : Structure containing the PhysarumSolver inputs
% * ListNodes             : The structure containing all the nodes
% 
%% Author: Aram Vroom (2016)
% Email:  aram.vroom@strath.ac.uk

disp('Initialising Physarum...')

% Solver Parameters:
InitialisedInputs = struct(...
                'LinearDilationCoefficient',          options.LinearDilationCoefficient,  ...   % Linear dilation coefficient 'm'
                'EvaporationCoefficient',             options.EvaporationCoefficient,  ...      % Evaporation coefficient 'rho'
                'GrowthFactor',                       options.GrowthFactorVal,  ...             % Growth factor 'GF'
                'NumberOfAgents',                     options.NumberOfAgents,  ...              % Number of virtual agents 'N_agents'
                'RamificationProbability',            options.RamificationProbability, ...      % Probability of ramification 'p_ram'
                'RamificationWeight',                 options.RamificationWeight,  ...          % Weight on ramification 'lambda'
                'MaximumRadiusRatio',                 options.MaximumRadiusRatio,  ...          % Maximum radius of the veins
                'MinimumRadiusRatio',                 options.MinimumRadiusRatio,  ...          % Minimum radius of the veins
                'StartingRadius',                     options.StartingRadius,  ...              % The starting radius of the veins
                'RamificationAmount',                 options.RamificationAmount,  ...          % The number of nodes initially generated for the ramification
                'PossibleDecisions',                  {options.Cities},  ...                    % The list of possible cities
                'MaxConsecutiveVis',                  options.MaxConsecutiveVis,  ...           % Maximum number of consecutive visits to each city
                'MaxVisits',                          options.MaxVisits,  ...                   % Maximum number of visits to each city
                'RootAttrib',                         options.RootAttrib,   ...                 % Attributes of the root
                'NodeCheckBoundaries',                options.NodeCheckBoundaries  , ...        % The values used by the MyCreatedNodeCheck file
                'Generations',                        options.Generations, ...                  % The number of generations
                'Viscosity',                          options.Viscosity, ...                    % The fluid viscosity "mu"
                'MinCommonNodesThres',                options.MinCommonNodesThres,  ...         % The minimum number of nodes two decision sequences should have in common for a restart to occur
                'IfZeroLength',                       options.IfZeroLength, ...                 % Value assigned to the length if it's zero (to prevent flux = inf)
                'MaxChildFindAttempts',               options.MaxChildFindAttempts, ...         % Max number of attempts that will be done to find additional children for a node
                'CostFunction',                       fitnessfcn,  ...                          % Reference to the fitness function file
                'MyNodeAttributes',                   options.MyNodeAttributes, ...             % Reference to the file containing the NodeAttributes class
                'MyAttributeCalcFile',                options.MyAttributeCalcFile, ...          % Reference to the file performing the attribute calculations
                'NodeIDCheckFile',                    options.MyNodeIDCheck, ...                % Reference to the file that checks the feasibility of a node's child using solely the unique ID
                'CreatedNodeCheckFile',               options.MyCreatedNodeCheck, ...           % Reference to the file that checks the feasibility of a node's child after the child's structure has been created
                'BestChainFile',                      options.MyBestChainFile, ...              % Reference to the file that determines the best chain
                'EndConditionsFile',                  options.MyEndConditionsFile, ...          % Reference toe the file that checks whether the end conditions have been reached
                'EndConditions',                      {options.EndConditions}, ...              % End conditions used in the options.MyEndConditionsFile file
                'AttributeIDIndex',                   options.AttributeIDIndex, ...             % Index of the attributes that determine the unique ID
                'AdditionalInputs',                   {options.AdditonalInputs}, ...            % Variable that can be used to store any additional information required in one of the user's files.
                'Sets',                               {sets}, ...                               % The sets of the possible attributes incorporated in the ID
                'RootName',                           char(options.RootName), ...               % The name of the root
                'MinPickProbability',                 options.MinPickProbability, ...           % The minimum probability for a feasible node to be picked before the algorithm changes its method of choosing a child
                'GenerateGraphPlot',                  options.GenerateGraphPlot, ...            % Indicator as to whether the algorithm should generate a graph plot animation                  
                'GraphPlotFileName',                  char(options.GraphPlotFileName), ...      % Name of the file that the graph plot animation will be saved as
                'GenerateTreePlot',                   options.GenerateTreePlot, ...             % Indicator as to whether the algorithm should generate a tree plot
                'SaveHistory',                        options.SaveHistory, ...                  % Indicator as to whether the algorithm should save the history of the radius of each vein and the path of each agent throughout the simulation
                'LowMem',                             options.LowMem ...                        % Indicator as to whether the algorithm should use the low-memory version of searchin for new nodes, where 1 is defined as "yes". Using the low-memory version is slower, but requires less memory. 
            );
        
      
% Error Checking:

% Display error if the vectors with number of possible decisions, max. number of consecutive
% resonance orbits & the max. number of visists is not equal
if (length(InitialisedInputs.MaxConsecutiveVis)~=length(InitialisedInputs.PossibleDecisions))
    error('Check size of PossibleDecisions, MaxConsecutiveVis and MaxVisits')
end

%Check if the sets input has the correct size
setnames = fieldnames(sets);
for i = 1:length(setnames)
    if (length(sets.(char(setnames(i))))~=length(options.Cities))
        error(strcat('Each field within the sets structure should contain a Cx1 cell array,', ...
            ' where C is the number of cities'));
    end
end

% Check if the size of the AttributeIDIndex option is correct
if (length(InitialisedInputs.AttributeIDIndex)~=length(fieldnames(sets)))
    error('Check size of AttributeIDIndex and sets')
end

% Add nodes that can be selected to the Inputs structure
InitialisedInputs.PossibleListNodes = PossChilds(options.Cities, sets);

% Create the list of nodes
ListNodes = CreateListNodes(InitialisedInputs);

% Clear workspace
clearvars -except ListNodes InitialisedInputs
disp('Initialisation Complete.')

end

