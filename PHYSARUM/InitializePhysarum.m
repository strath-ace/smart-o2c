function [InitalizedInputs, ListNodes] = InitializePhysarum(UserInputs)
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

%Retrieve the user settings
PhysarumSettings

%Solver Parameters:
InitalizedInputs = struct('LowThrust',                UserInputs.LowThrust,  ... %Set to 1 for low-thrust, 0 for high-thrust
                'LinearDilationCoefficient',          UserInputs.LinearDilationCoefficient,  ... %Linear dilation coefficient 'm'
                'EvaporationCoefficient',             UserInputs.EvaporationCoefficient,  ... %Evaporation coefficient 'rho'
                'GrowthFactor',                       UserInputs.GrowthFactor,  ... %Growth factor 'GF'
                'NumberOfAgents',                     UserInputs.NumberOfAgents,  ... %Number of virtual agents 'N_agents'
                'RamificationProbability',            UserInputs.RamificationProbability, ... %Probability of ramification 'p_ram'
                'RamificationWeight',                 UserInputs.RamificationWeight,  ... %Weight on ramification 'lambda'
                'MaximumRadiusRatio',                 UserInputs.MaximumRadiusRatio,  ... %Maximum radius of the veins
                'MinimumRadiusRatio',                 UserInputs.MinimumRadiusRatio,  ... %Minimum radius of the veins
                'StartingRadius',                     UserInputs.StartingRadius,  ... %The starting radius of the veins
                'RamificationAmount',                 UserInputs.RamificationAmount,  ... %The number of nodes initially generated for the ramification
                'PossibleDecisions',                  {PossibleDecisions},  ... %The list of possible targets
                'MaxConsecutiveRes',                  {MaxConsecutiveRes},  ... %Maximum number of consecutive resonance orbits to each target
                'MaxVisits',                          {MaxVisits},  ... %Maximum number of visits to each target
                'RootChar',                           UserInputs.RootChar,   ... %Characteristic of the root
                'Generations',                        UserInputs.Generations, ... %The number of generations
                'Viscosity',                          UserInputs.Viscosity, ... %The viscocity of the "fluid" 
                'DeterminingCharacteristic',          UserInputs.DeterminingCharacteristic, ... %The index of the determining characteristic in the 'characteristics' field
                'MinCommonNodesThres',                UserInputs.MinCommonNodesThres  ... %The minimum number of nodes two decision sequences should have in common for a restart to occur
                );               

%Create a list of all possible characteristics & decisions
possiblecharacteristics = mincharboundary:stepsize:maxcharboundary;

%Create list of nodes that can be selected
for i = 1:(length(possiblecharacteristics))
    for j = 1:length(PossibleDecisions)
        possdeccharvec(i,j) = strcat(PossibleDecisions(j),'_',num2str(possiblecharacteristics(i)));
    end
end
possdeccharvec = reshape(possdeccharvec, [1, numel(possdeccharvec)]);

%Add nodes that can be selected to the Inputs structure
InitalizedInputs.PossibleListNodes = possdeccharvec;


%Display error if the vectors with number of possible decisions, max. number of consecutive
%resonance orbits & the max. number of visists is not equal
if ~(max(size(MaxConsecutiveRes))==max(size(PossibleDecisions)))
    error('Check size of PossibleDecisions, MaxConsecutiveRes and MaxVisits')
end

%Create the list of nodes
ListNodes = CreateListNodes(InitalizedInputs);

end

