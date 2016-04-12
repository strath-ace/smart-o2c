% This script contains the settings for the physarum solver. As such, this
% script will be run before the solver itself is started
%
% Inputs:
% * 
%
% Outputs: 
% * UserInputs         : Structure containing the uesr's PhysarumSolver inputs
% * PossibleDecisions  : The possible decisions (targets) that the solver
%                        can choose from
% * MaxConsecutiveRes  : The maximum number of resonance orbits to each
%                        target (set to -1 to ignore)
% * MaxVisits          : The maximum nubmer of visists to each target (set
%                        to -1 to ignore)
% * mincharboundary    : The minimum boundary for the characteristic(s)
% * maxcharboundary    : The maximum boundary for the characteristic(s)
% * stepsize           : The stepsize with which the characteristic(s) is/are 
%                        to be evaluated
%
% Author: Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk

PossibleDecisions = {'A', 'B', 'C', 'D', 'E', 'F'}; %The possible decisions (targets) that the solver can choose from
MaxConsecutiveRes = -1*ones(1, 6); %The maximum number of resonance orbits to each target (set to -1 to ignore)
MaxVisits = [3 1 1 1 1 1]; %The maximum nubmer of visists to each target (set to -1 to ignore)

%Set boundaries. These are used to generate a list of all the possible
%nodes that can be chosen from.
mincharboundary  = 0; %The minimum boundary of the characteristic
maxcharboundary = 100; %The maximum boundary of the characteristic
stepsize = 1; %The step size by which the characteristic is evaluated


UserInputs = struct('LowThrust',                          1,  ... %Set to 1 for low-thrust, 0 for high-thrust
                    'LinearDilationCoefficient',          20,  ... %Linear dilation coefficient 'm'
                    'EvaporationCoefficient',             0,  ... %Evaporation coefficient 'rho'
                    'GrowthFactor',                       0,  ... %Growth factor 'GF'
                    'NumberOfAgents',                     3,  ... %Number of virtual agents 'N_agents'
                    'RamificationProbability',            0.1, ... %Probability of ramification 'p_ram'
                    'RamificationWeight',                 1,  ... %Weight on ramification 'lambda'
                    'MaximumRadius',                      50,  ... %Maximum radius of the veins
                    'MinimumRadius',                      1e-3,  ... %Minimum radius of the veins
                    'StartingRadius',                     1,  ... %The starting radius of the veins
                    'RamificationAmount',                 5,  ... %The number of nodes initially generated for the ramification
                    'RootChar',                           0,  ... %Characteristic of the root
                    'Generations',                        3  ... %The number of generations
                );  

                