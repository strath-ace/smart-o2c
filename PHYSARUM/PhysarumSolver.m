function [Inputs] = PhysarumSolver()
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


%Settings as to which targets are available, how many max. consecutive resonance
%orbits one wants and the maximum number of visits to each target
PossibleDecisions = {'A','B','C','D','E','F'};
MaxConsecutiveRes = -1*ones(1,6);
MaxVisits = [3 1 1 1 1 1];

if ~(max(size(MaxConsecutiveRes))==max(size(PossibleDecisions)))
    error('Check size of PossibleDecisions, MaxConsecutiveRes and MaxVisits')
end

%Set boundaries
mincharboundary  = 0;
maxcharboundary = 9;
stepsize = 1;

%Solver Parameters:
Inputs = struct(...
'LowThrust',                          1,  ... %Set to 1 for low-thrust, 0 for high-thrust
'LinearDilationCoefficient',          0,  ... %Linear dilation coefficient 'm'
'EvaporationCoefficient',             0,  ... %Evaporation coefficient 'rho'
'GrowthFactor',                       0,  ... %Growth factor 'GF'
'NumberOfAgents',                     0,  ... %Number of virtual agents 'N_agents'
'RamificationProbability',            0.2, ... %Probability of ramification 'p_ram'
'RamificationWeight',                 0,  ... %Weight on ramification 'lambda'
'MaximumRadius',                      0,  ... %Maximum radius of the veins
'MinimumRadius',                      0,  ... %Minimum radius of the veins
'StartingRadius',                     0,  ... %The starting radius of the veins
'RamificationAmount',                 5,  ... %The number of nodes initially generated for the ramification
'PossibleDecisions',{PossibleDecisions},  ... %The list of possible targets
'MaxConsecutiveRes',{MaxConsecutiveRes},  ... %Maximum number of consecutive resonance orbits to each target
'MaxVisits',                {MaxVisits},  ... %Maximum number of visits to each target
'RootChar',                           0   ... %Characteristic of the root
);               

%Create a list of all possible characteristics & decisions
possiblecharacteristics = mincharboundary:stepsize:maxcharboundary;

%Create list of nodes that can be selected
for i = 1:(length(possiblecharacteristics))
    for j = 1:length(PossibleDecisions)
    possdeccharvec(i,j) = strcat(PossibleDecisions(j),'_',num2str(possiblecharacteristics(i)));
    end
end
possdeccharvec = reshape(possdeccharvec,[1,numel(possdeccharvec)]);

Inputs.PossibleNodes = possdeccharvec;
end

