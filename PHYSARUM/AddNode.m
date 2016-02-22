function [Nodes] = AddNode(Inputs,Nodes,parent,linkstoadd)
% This function adds a node to the existing graph
%   
% Inputs:
% * Inputs : Structure containing the PhysarumSolver inputs
% * Nodes  : Structure that contains the currently existing nodes
% * parent : The parent of the node to be added
% * linkstoadd: vector containing the rows of the nodes that the node to be 
%   added is connected to
%
% OUTPUTS: 
% * Nodes: Previous Nodes structure but with additional node
%
% Author: Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk

%Find the row to which elements should be added
addrow = length(Nodes.identifier)+1;

%Define a variable that contains the number of links that are added, for
%ease of writing
numberoflinks = length(linkstoadd);

%Add the (starting) characteristics of each vein to the Nodes structure
Nodes.identifier{addrow,1}        = addrow;
Nodes.parent{addrow,1}            = parent;
Nodes.links{addrow,1}             = linkstoadd;
Nodes.radius{addrow,1}            = ones(1,numberoflinks)*Inputs.StartingRadius;
Nodes.pressure_gradient{addrow,1} = ones(1,numberoflinks);
Nodes.lengths{addrow,1}           = ones(1,numberoflinks);
Nodes.fluxes{addrow,1}            = zeros(1,numberoflinks);
Nodes.probabilities{addrow,1}     = ones(1,numberoflinks);

 