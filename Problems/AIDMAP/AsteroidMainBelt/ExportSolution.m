function [] = ExportSolution(InitializedInputs,ListNodes,Nodes,filename)
% This script exports the found solution to a file
%
% Inputs:
% * InitializedInputs  : The structure containing the options set by the
%                        user
% * ListNodes          : Structure containing the initial list of nodes
% * Nodes              : The nodes to be saved as a solution
%
% Outputs: 
% *
%
% Author: Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk

%Loop over the nodes
for i = 1:length(Nodes)
    
    %Save each node as a field in the SolutionNodes structure
    SolutionNodes.(char(Nodes(i))) = ListNodes.(char(Nodes(i)));
end

%Save the initialized inputs and the solution
save(filename,'SolutionNodes','InitializedInputs');
end

