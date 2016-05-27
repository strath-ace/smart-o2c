function [x,fval,exitflag,output] = optimisation_aidmap(fitnessfcn,sets,options)
%function [output] = optimisation_aidmap(fitnessfcn,sets,options)

%% optimisation_aidmap: This function handles the optimisation using the AIDMAP algorithm
% Extensive function description
% (If you need to insert formulas use latex conventions: 
% $x_1+x_2$) 
%% Inputs:
%
% * fitnessfcn : The function reference to the fitness function
% * sets       : The sets of the possible attributes
% * options    : Structure containing the options sets by the user
%
%
%% Output:
% * x      : The best chain of nodes found
% * fval   : The cost of the best chain of nodes found
% * output : The structure containing the AIDMAP algorithm's outputs
%
% Author: Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk

tic;

[InitializedInputs,ListNodes] = InitializePhysarum(fitnessfcn,options,sets);

[output.Solutions, BestSolution, InitializedInputs, output.ListNodes, output.Agents] = PhysarumSolver(InitializedInputs, ListNodes);

x = BestSolution.BestChain;
fval = BestSolution.BestCost;

output.CompTime = toc;

exitflag = 0;

end
