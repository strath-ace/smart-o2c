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
MaxConsecutiveRes = [0 0 0 0 0 0];
MaxVisits = [3 1 1 1 1 1];

if ~(max(size(MaxConsecutiveRes))==max(size(PossibleDecisions)))
    error('Check size of PossibleDecisions, MaxConsecutiveRes and MaxVisits')
end


%Solver Parameters:
Inputs = struct(...
'LowThrust',                          1,  ... %Set to 1 for low-thrust, 0 for high-thrust
'LinearDilationCoefficient',          0,  ... %Linear dilation coefficient 'm'
'EvaporationCoefficient',             0,  ... %Evaporation coefficient 'rho'
'GrowthFactor',                       0,  ... %Growth factor 'GF'
'NumberOfAgents',                     0,  ... %Number of virtual agents 'N_agents'
'RamificationProbability',            0,  ... %Probability of ramification 'p_ram'
'RamificationWeight',                 0,  ... %Weight on ramification 'lambda'
'MaximumRadius',                      0,  ... %Maximum radius of the veins
'MinimumRadius',                      0,  ... %Minimum radius of the veins
'StartingRadius',                     0,  ... %The starting radius of the veins
'RamificationAmount',                 5,  ... %The number of nodes initially generated for the ramification
'PossibleDecisions',{PossibleDecisions},  ...       
'MaxConsecutiveRes',{MaxConsecutiveRes},  ...
'MaxVisits',                {MaxVisits}  ...
);               

end

