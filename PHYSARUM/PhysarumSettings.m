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


%Targets the Physarum can choose from
Targets = {'neo163693', 'neo164294', 'neo1998DK36', 'neo2004JG6', 'neo2005TG45','neo2006WE4','neo2007EB26','neo2008EA32' ,'neo2008UL90' ,'neo2010XB11' ,'neo2012VE46' ,'neo2013JX28'}; 

%The maximum number of resonance orbits to each target (set to -1 to ignore)
MaxConsecutiveRes = -1*ones(1, length(Targets)); 

%The maximum nubmer of visists to each target (set to -1 to ignore)
MaxVisits = ones(1, length(Targets)); %[3 1 1 1 1]; 

%Index of the attributes that determine the unique ID (characteristics)
AttributeIDIndex = [1 3];  

%Set values of these characteristics. If the values are the same for every
%target, a vector can be used. Otherwise a matrix should be set with the
%characteristics for every target, where every row denotes a target. If the
%number of values a characteristic can take differs per target (thus making 
%it impossible to define a matrix), a cell array can be used for that characteristic 
%as well.
charvalues{1} = 0:1:20;
charvalues{2} = {{2 3 4}; {2 3}};

%Low Thrust Flag. Set to 1 for low-thrust, 0 for high-thrust
LowThrust = 1;

%Linear dilation coefficient 'm'
LinearDilationCoefficient = 20;

%Evaporation coefficient 'rho'
EvaporationCoefficient = 0;

%Growth factor 'GF'
GrowthFactorVal = 0;

%Number of virtual agents 'N_agents'
NumberOfAgents = 3;

%Probability of ramification 'p_ram'
RamificationProbability = 0.15;

%Weight on ramification 'lambda'
RamificationWeight = 1;

%Maximum ratio between the link's radius & the starting radius
MaximumRadiusRatio = 1000;

%Maximum ratio between the link's radius & the starting radius
MinimumRadiusRatio = 1e-3;

%The starting radius of the veins
StartingRadius = 1;

%The number of nodes initially generated for the ramification
RamificationAmount = 5;   

%Characteristic of the root
RootChar = [0 0];   

%The number of generations
Generations = 5;   

%The viscocity of the "fluid" 
Viscosity = 1;  

%The index of the determining characteristic in the 'characteristics' field
DeterminingCharacteristic = 1;  

%The minimum number of nodes two decision sequences should have in common for a restart to occur
MinCommonNodesThres = 7;  

%Value assigned to the length if it's zero (to prevent flux = inf)
IfZeroLength = 1e-15; 

%The functio nreference to the cost function
CostFunction = @MyCostFunction; 

%The class that contains the node attributes
NodeAttributes = @MyAttributes; 

%The project directory
ProjectDirectory = 'C:\Users\ckb16114\Desktop\Internship\Code\Developing\Atira Algorithm';





                
