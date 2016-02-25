function [newNode] = CreateNode(node_ID,parent)
%This function creates the structure for a new node.
%
% Inputs:
% * node_ID : a string containing the node_ID of the node to be generated
% * parent  : a string containing the node_ID of the parent node
%
% Outputs: 
% * newNode : a structure describing the new node and its attributes
%
% Author: Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk

newNode = struct(node_ID, struct(                        ...
                     'parent',            parent,...
                     'children',          [],... % Matrix that holds the nodes' connections to each other
                     'radius',            [],... % The radius of each connection
                     'pressure_gradient', [],... % The pressure gradient over each connection
                     'lengths',           [],... % The length of each connection
                     'fluxes',            [],... % Matrix containing each connection's flux
                     'probabilities',     [] ... % Matrix containing the probability for each connection
                     ));

end

