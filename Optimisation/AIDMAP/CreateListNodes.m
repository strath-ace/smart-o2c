function [ListNodes] = CreateListNodes(Inputs)
% This function creates the structure describing the graph.
%
% Inputs:
% * Inputs : Structure containing the PhysarumSolver inputs
%
% Outputs: 
% * ListNodes  : Empty structure that will contain the Graph's information
%
% Author: Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk

%Create the ListNodes structure
ListNodes = struct(struct('Root',   ...
                         struct(...
                         'node_ID',           'Root',...
                         'parent',            [],... % The parent of the node
                         'children',          [],... % Matrix that holds the nodes' connections to each other
                         'radius',            [],... % The radius of each connection
                         'length',           [],... % The length of each connection
                         'flux',            [],... % Matrix containing each connection's flux
                         'attributes',   [SetNodeAttributes(Inputs,Inputs.RootAttrib)],... % Attributes that describe this node (such as orbital elements & ToF .)
                         'previousdecisions', [],... % List of the previous decisions made
                         'possibledecisions', {Inputs.PossibleDecisions}, ... % Targets that can still be visisted by the node
                         'VisitsLeft',        {Inputs.MaxVisits} ... % Vector containing the number of times each target cna still be visisted
                     )));
           
end
