function [x, fval, exitflag, output] = optimise_aidmap(fitnessfcn, sets, options)
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%-----------Copyright (C) 2016 University of Strathclyde-------------
%
%
%
%% optimisation_aidmap:
% The AIDMAP algorithm is a combinatorial algorithm that takes its inspiration 
% from the Physarum Polycephalum mould. To simulate this mould, a number of virtual agents 
% are used to resemble nutrients inside veins that move through and allow 
% the incremental branching of new veins. Each of these veins have a radius, 
% length and an amount of flux going through them, where the latter depends on the first two. 
% These characteristics are stored in so-called nodes that are placed in between every two veins. 
% Aside from the vein characteristics, each node also contains the problem-specific attributes. 
% As each node resembles a certain decision, a set of nodes connected by veins can resemble
% the set of consecutive decisions present in discrete decision making problems.
%
%% Inputs:
% * fitnessfcn : Function handle to cost function (real function)
% * sets       : The structure containing the values that the optimisation
%                variables can have for each city (where "city" is defined as done 
%                in the traveling salesman problem). Thus, each field
%                contains a Cx1 cell array, where C is the number of cities and
%                each cell in turn contains the discrete set of values the optimisation
%                variable can have
% * options    : Structure containing the options sets by the user. 
%                Note that all these options need to be defined. 
%                The references to files starting with "My" are files that
%                are to be provided by the user. Examples of these files
%                and their structure can be found in the Problems folder of
%                the SMART-O2C Toolbox
%               * options.LinearDilationCoefficient: Linear dilation
%                   coefficient 'm' [real number]
%               * options.EvaporationCoefficient: Evaporation coefficient
%                   'rho' [real number]
%               * options.GrowthFactorVal: Growth factor 'GF' [real number]
%               * options.NumberOfAgents: Number of virtual agents
%                   'N_agents' [integer]
%               * options.RamificationProbability: Probability of
%                   ramification 'p_ram' [real number between 0 and 1, where
%                   1 is a 100% probability for an agent to ramificate]
%               * options.RamificationWeight: Weight on ramification
%                   'lambda' [real number, where a larger value puts more
%                   weight on ramification]
%               * options.MaximumRadiusRatio: Maximum radius of the veins [real number]
%               * options.MinimumRadiusRatio: Minimum radius of the veins [real number]
%               * options.StartingRadius: The starting radius of the veins [real number]
%               * options.RamificationAmount: The number of nodes initially
%                   generated for the ramification [integer]
%               * options.Cities: The list of possible cities [1xC string
%                   array, where C is the number of cities]
%               * options.MaxConsecutiveVis: Maximum number of consecutive
%                   visits to each city. Set maxima to -1 if no maximum
%                   defined [1xC vector of integers, with C being the number of cities]
%               * options.MaxVisits: Maximum number of visits to each
%                   city. Set maxima to -1 if no maximum defined
%                   [1xC vector of integers, with C being the number of cities]
%               * options.RootAttrib: Values of the optimisation variables
%                   at the root [1xV vector of real numbers, where V is the
%                   number of optimisation variables]
%               * options.NodeCheckBoundaries: The values of the boundaries that
%                   define whether a link created to a new node is valid. Used
%                   in the user's MyCreateNodeCheck file [vector of real
%                   numbers]
%               * options.Generations: The number of generations [integer]
%               * options.Viscosity: The fluid viscosity "mu" [real number]
%               * options.MinCommonNodesThres: The minimum number of nodes
%                   two agents in a generation should have in common for a
%                   restart to occur [integer]
%               * options.IfZeroLength: Value assigned to the length if
%                   it's zero (to prevent flux = inf) [real number]
%               * options.MaxChildFindAttempts: Max number of attempts that
%                   will be done to find additional children for a node
%                   [integer]
%               * options.MyNodeAttributes:  Reference to the file
%                   containing the NodeAttributes class, which defines the
%                   problems-specific attributes each node has          
%               * options.MyAttributeCalcFile: Reference to the file that
%                   calcules the problem-specific attributes
%               * options.MyNodeIDCheck: Reference to the file that checks
%                   the feasibility of a node's child using solely the
%                   optimisation variables. Allows the algorithm to determine
%                   whether a connection to a child node is valid, before
%                   generating the full structure of the child node
%               * options.MyCreatedNodeCheck: Reference to the file that checks 
%                   the feasibility of a node's child after the child's 
%                   structure has been created
%               * options.MyBestChainFile: Reference to the file that
%                   determines the best chain (decision path)
%               * options.MyEndConditionsFile: Reference to the file that checks 
%                   whether the end conditions have been reached
%               * options.EndConditions: End conditions used in the
%                   options.MyEndConditionsFile file [cell array]
%               * options.AttributeIDIndex: Index of the optimisation variables in the 
%                   MyNodeAttributes class [1xV vector of integer, 
%                   where V is the number of optimisation variables]
%               * options.AdditonalInputs: Variable that can be used to
%                   store any additional information required in one of the
%                   user's files. [cell array]
%               * options.RootName: The name of the root [string]
%               * options.MinPickProbability: The minimum probability for a
%                   feasible node to be picked from the list of all possible children
%                   before the algorithm changes its method of choosing a
%                   child [real number between 0 and 1, where setting it to 1
%                   makes the algorithm immediately use the alternative
%                   method of picking a child (this alternative method is
%                   slower, but is more capable of finding valid childs
%                   when the probability of picking a valid child from the
%                   full list of possible childs is very low)
%               * options.GenerateGraphPlot: Indicator as to whether the
%                   algorithm should generate a graph plot animation, where 1 is
%                   defined as 'yes'
%               * options.GraphPlotFileName: Name of the file that the
%                   graph plot animation will be saved as [string]
%               * options.GenerateTreePlot: Indicator as to whether the
%                   algorithm should generate a tree plot, where 1 is defined
%                   as "yes"
%               * options.SaveHistory: Indicator as to whether the
%                   algorithm should save the history of the radius of each
%                   vein and the path of each agent throughout the simulation, 
%                   where 1 is defined as 'yes'
%               * options.LowMem = 0: Indicator as to whether the algorithm 
%                   should use the low-memory version of searching for new nodes, 
%                   where 1 is defined as "yes". Using the low-memory version is slower
%                   for small problems, but requires less memory. For large cases, the LowMem
%                   version may run faster due to the fact that certain node-selection
%                   variables are not not saved when this version is used
%
%% Outputs:
% * x           : The best chain of nodes found [string array of node IDs]
% * fval        : The cost of the best chain of nodes found [real value]
% * exitflag    : Indicator as to whether a solution has been found. 
%                 Is 1 if a solution has been found, and 0 otherwise
% * output      : The structure containing the AIDMAP algorithm's outputs
%                * output.Solutions : The structure containing the solutions found
%                   * output.Solutions.Nodes: cell array containing all the
%                      solutions (paths) found
%                   * output.Solutions.Costs: cell array containing the costs
%                      corresponding to each link of each solution found
%                * output.ListNodes: Structure containing the final structure with the nodes
%                * output.Agents: the structure containing the set of agents and their characteristics
%                * output.History: the vein radii and agent movement
%                   throughout the generations
%                * output.funcalls: number of cost function calls
%                * output.CompTime: the computation time of the algorithm
%
%% Author(s): Aram Vroom (2016), Juan Manuel Romero Martin (2014) and Luca Masi (2012)
% Email: aram.vroom@strath.ac.uk juan.romero-martin@strath.ac.uk luca.masi@strath.ac.uk
%
%% References:
% * Masi, Luca and Vasile, Massimiliano (2014) A multidirectional Physarum solver for the automated design of space trajectories. In: Proceedings of the 2014 IEEE Congress on Evolutionary Computation, CEC 2014. Institute of Electrical and Electronics Engineers Inc., pp. 2992-2999. ISBN 978-1-4799-6626-4
% https://pure.strath.ac.uk/portal/files/39689086/Masi_Vasile_CEC_2014_A_multidirectional_physarum_solver_automated_design_space_trajectories_Jul_2014.html


% Start clock
tic;

% Initialise the AIDMAP algorithm
[InitialisedInputs, ListNodes] = InitialisePhysarum(fitnessfcn, options, sets);

% Run the algorithm
[output.Solutions, BestSolution, output.ListNodes, output.Agents, output.History, output.funccalls] = PhysarumSolver(InitialisedInputs, ListNodes);

% Retrieve the best solution
x = BestSolution.BestChain;
fval = BestSolution.BestCost;

% Exitflag = 0 denotes no solution found
exitflag = 0;

% If the length of x is more than 1, a solution was found
if (length(x{1})>1)
    
    % Set the exitflag to 1
    exitflag = 1;
    
    % If the user has requested the generation of the graph plot animation, 
    % generate this plot
    if (options.GenerateGraphPlot==1)
        PhysarumGraphPlot(options, output.ListNodes, output.History);
    end
    
    % If the user has requested the generation of the tree plot, generate this plot
    if (options.GenerateTreePlot==1)
        PhysarumTreePlot(output.ListNodes)
    end
end

% End clock
output.CompTime = toc;

end
