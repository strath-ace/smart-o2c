function [Nodes] = CreateGraph(Inputs)
% This function creates the structure describing the graph.
%
% Inputs:
% * none
%
% Outputs: 
% * Nodes : Empty structure that will contain the Graph's information
%
% Author: Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk

%Create the Nodes structure
Nodes = struct('Root',struct(                    ...
                     'node_ID',           'Root',...
                     'parent',            [],... % The parent of the node
                     'children',          [],... % Matrix that holds the nodes' connections to each other
                     'radius',            [],... % The radius of each connection
                     'pressure_gradient', [],... % The pressure gradient over each connection
                     'lengths',           [],... % The length of each connection
                     'fluxes',            [],... % Matrix containing each connection's flux
                     'probabilities',     [],... % Matrix containing the probability for each connection
                     'characteristics',   [],... % Characteristics that describe this node (such as orbital elements & ToF .)
                     'previousdecisions', [],... % List of the previous decisions made
                     'possibledecisions', [],... % Targets that can still be visisted by the node
                     'VisitsLeft',        [Inputs.MaxVisits] ... % Vector containing the number of times each target cna still be visisted
                 ),                         ...
               'ListNodes',struct('Root',struct(...
                     'node_ID',           'Root',...
                     'parent',            [],... % The parent of the node
                     'children',          [],... % Matrix that holds the nodes' connections to each other
                     'radius',            [],... % The radius of each connection
                     'pressure_gradient', [],... % The pressure gradient over each connection
                     'lengths',           [],... % The length of each connection
                     'fluxes',            [],... % Matrix containing each connection's flux
                     'probabilities',     [],... % Matrix containing the probability for each connection
                     'characteristics',   [Inputs.RootChar],... % Characteristics that describe this node (such as orbital elements & ToF .)
                     'previousdecisions', [],... % List of the previous decisions made
                     'possibledecisions', {Inputs.PossibleDecisions}, ... % Targets that can still be visisted by the node
                     'VisitsLeft',        {Inputs.MaxVisits} ... % Vector containing the number of times each target cna still be visisted
                 )));
           
end
