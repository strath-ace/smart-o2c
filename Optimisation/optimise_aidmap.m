function [x,fval,exitflag,output] = optimise_aidmap(fitnessfcn,sets,options)

%% optimisation_aidmap: 
% the AIDMAP algorithm is a combinatorial algorithm that takes its inspiration from the Physarum Polycephalum mould. To simulate this mould, a number of virtual agents are used to resemble nutrients inside veins that move through and allow the incremental branching of new veins. Each of these veins have a radius, length and an amount of flux going through them, where the latter depends on the first two. These characteristics are stored in so-called nodes that are placed in between every two veins. Aside from the vein characteristics, each node also contains the problem-specific attributes. As each node resembles a certain decision, a set of nodes connected by veins can resemble the set of consecutive decisions present in the aforementioned discrete decision making problems.
%
%% Inputs:
%
% * fitnessfcn : function handle to cost function (real function)
% * sets       : set of nodes bounds (is a matrix dx2 where d is the dimension of the discrete variables asscoiated to each node)
% * options    : Structure containing the options sets by the user
%               * options.LinearDilationCoefficient: Linear dilation coefficient 'm'
%               * options.EvaporationCoefficient: Evaporation coefficient 'rho'
%               * options.GrowthFactorVal: Growth factor 'GF'
%               * options.NumberOfAgents: Number of virtual agents 'N_agents'
%               * options.RamificationProbability: Probability of ramification 'p_ram'
%               * options.RamificationWeight: Weight on ramification 'lambda'
%               * options.MaximumRadiusRatio: Maximum radius of the veins
%               * options.MinimumRadiusRatio: Minimum radius of the veins
%               * options.StartingRadius: The starting radius of the veins
%               * options.RamificationAmount: The number of nodes initially generated for the ramification
%               * {options.Targets}: The list of possible targets
%               * options.MaxConsecutiveRes: Maximum number of consecutive resonance orbits to each target
%               * options.MaxVisits: Maximum number of visits to each target
%               * options.RootAttrib: Attributes of the root
%               * options.NodeCheckBoundaries: The values used by the MyCreatedNodeCheck file
%               * options.Generations: The number of generations
%               * options.Viscosity: The viscocity of the "fluid" 
%               * options.MinCommonNodesThres: The minimum number of nodes two decision sequences should have in common for a restart to occur
%               * options.IfZeroLength: Value assigned to the length if it's zero (to prevent flux = inf)
%               * options.MaxChildFindAttempts: Max number of attempts that will be done to find additional children for a node
%               * options.NodeAttributes: Reference to the file containing the NodeAttributes class
%               * options.MyAttributeCalcFile: Reference to the file performing the attribute calculations
%               * options.MyNodeIDCheck: Reference to the file that checks the feasibility of a node's child using solely the unique ID
%               * options.MyCreatedNodeCheck: Reference to the file that checks the feasibility of a node's child after the child's structure has been created
%               * options.MyBestChainFile: Reference to the file that determines the best chain
%               * options.MyEndConditionsFile: Reference toe the file that checks whether the end conditions have been reached
%               * {options.EndConditions}: End conditions used in the options.MyEndConditionsFile file
%               * options.AttributeIDIndex: Index of the attributes that determine the unique ID
%               * {options.AdditonalInputs}: Made into cell in case multiple additional outpu
%               * options.RootName: The name of the root
%               * options.MinPickProbability: The minimum probability for a feasible node to be picked before the algorithm changes its method of choosing a child
%               * options.GenerateGraphPlot: 
%               * options.SaveHistory: 
%
%% Output:
% * x           : The best chain of nodes found (best integer solution)
% * fval        : The cost of the best chain of nodes found
% * exitflag    : =1 if solution has been found, =0 otherwise
% * output      : The structure containing the AIDMAP algorithm's outputs
%               * output.Solutions: The structure containing the solutions found
%               * output.ListNodes: Structure containing the final structure with the nodes
%               * output.Agents: the structure containing the set of agents and their characteristics
%               * output.History: 
%               * output.funcalls: number of cost function calls
%               * output.CompTime: 
%
%% Author(s): Aram Vroom (2016)
% Email:  aram.vroom@strath.ac.uk
%
%% References


%Start clock
tic;

%Initialize the AIDMAP algorithm
[InitializedInputs,ListNodes] = InitializePhysarum(fitnessfcn,options,sets);

%Run the algorithm
[output.Solutions, BestSolution, InitializedInputs, output.ListNodes, output.Agents, output.History, output.funccalls] = PhysarumSolver(InitializedInputs, ListNodes);

%Retrieve the best solution
x = BestSolution.BestChain(1);
fval = BestSolution.BestCost(1);

%Save the options
output.options = InitializedInputs;

%End clock
output.CompTime = toc;

exitflag = 1;
if (length(x{1})==1)
    exitflag = 0; %Exitflag = 0 denotes no solution found
end

end
