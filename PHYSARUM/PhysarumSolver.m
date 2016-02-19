function [ output_args ] = PhysarumSolver()
% This is the main file of the physarum solver. 
% It contains the logic for the algorithm and the solver parameters.
%
% Inputs:
% * 
%
% Outputs: 
% * 
%
% Author: Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk


%Solver Parameters:
Inputs = struct(...
'LowThrust',                    1,  ... %Set to 1 for low-thrust, 0 for high-thrust
'LinearDilationCoefficient',    0,  ... %Linear dilation coefficient 'm'
'EvaporationCoefficient',       0,  ... %Evaporation coefficient 'rho'
'GrowthFactor',                 0,  ... %Growth factor 'GF'
'NumberOfAgents',               0,  ... %Number of virtual agents 'N_agents'
'RamificationProbability',      0,  ... %Probability of ramification 'p_ram'
'RamificationWeight',           0   ... %Weight on ramification 'lambda'
);


            
            

end

