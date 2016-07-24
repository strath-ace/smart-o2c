function [] = ExportSolution(ListNodes, Nodes, filename)
%% ExportSolution: This script exports the found solution to a file
% 
%% Inputs:
% * ListNodes          : Structure containing the initial list of nodes
% * Nodes              : The nodes to be saved as a solution
% * filename           : The name of the file the solution is saved as
% 
%% Outputs: 
% *
% 
%% Author(s): Aram Vroom (2016)
% Email:  aram.vroom@strath.ac.uk

% Loop over the nodes
for i = 1:length(Nodes)
    
    % Save each node as a field in the SolutionNodes structure
    SolutionNodes.(char(Nodes(i))) = ListNodes.(char(Nodes(i)));
end

% Save the solution
save(strcat(filename, '.mat'), 'SolutionNodes');
end

