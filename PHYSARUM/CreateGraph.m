function [Nodes] = CreateGraph()
%This function creates the structure describing the graph
%
% INPUTS:
%   none
%
% OUTPUTS: 
%   Nodes      - Empty structure that will contain the Graph's / Tree's
%                information
%
% CREATED BY:
%   Aram Vroom - 2016

Nodes = struct('identifier',      [1],... %
               'parent',          [0],... %Matrix containing each node's parent
               'links',           [0],... %Matrix that holds the nodes' connections to each other
               'costs',           [0],... %Matrix containing each connection's cost
               'probabilities',   [0] ... %Matrix containing the probability for each connection
               ); 
 end

