function [] = ExportSolution(ListNodes, Nodes, filename)
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
<<<<<<< HEAD
%-----------Copyright (C) 2018 University of Strathclyde and Authors-----------
=======
%-----------Copyright (C) 2016 University of Strathclyde-------------
>>>>>>> 5b7361d93c9119cf1d2e9e6c885bed93f924d71b
%
%
%
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

