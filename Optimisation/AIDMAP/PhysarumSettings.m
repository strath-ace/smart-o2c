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
% * minattribboundary    : The minimum boundary for the attribute(s)
% * maxattribboundary    : The maximum boundary for the attribute(s)
% * stepsize           : The stepsize with which the attribute(s) is/are 
%                        to be evaluated
%
% Author: Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk


%Targets the Physarum can choose from
Targets = {'neo163693', 'neo164294', 'neo1998DK36', 'neo2004JG6', 'neo2005TG45','neo2006WE4','neo2007EB26','neo2008EA32' ,'neo2008UL90' ,'neo2010XB11' ,'neo2012VE46' ,'neo2013JX28'}; 
%targets = {'A','B','C'};

%The maximum number of resonance orbits to each target (set to -1 to ignore)
MaxConsecutiveRes = 0*ones(1, length(Targets)); 

%The maximum nubmer of visists to each target (set to -1 to ignore)
MaxVisits = 1*ones(1, length(Targets)); %[3 1 1 1 1]; 

%Index of the attributes that determine the unique ID
AttributeIDIndex = [11 10];  

%Set values of these attributes. If the values are the same for every
%target, a vector can be used. Otherwise a matrix should be set with the
%attributes for every target, where every row denotes a target. If the
%number of values an attribute can take differs per target (thus making 
%it impossible to define a matrix), a cell array can be used for that attribute 
%as well.
attribvalues{1} = 30:10:365;
load('epochsnode.mat')
attribvalues{2} = epochsnode;

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

%Attributes of the root
RootAttrib = [0 0];   

%The number of generations
Generations = 5;   

%The viscocity of the "fluid" 
Viscosity = 1;  

%The index of the determining attribute in the 'attributes' field
DeterminingAttribute = 1;  

%The minimum number of nodes two decision sequences should have in common for a restart to occur
MinCommonNodesThres = 7;  

%Value assigned to the length if it's zero (to prevent flux = inf)
IfZeroLength = 1e-15; 

%The function reference to the cost function
CostFunction = @MyCostFunction; 

%The class that contains the node attributes
NodeAttributes = @MyAttributes; 

%The project directory
ProjectDirectory = 'C:\Users\ckb16114\Desktop\Internship\Code\Developing\Atira Algorithm';





                
