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

Nodes = struct('identifier',      [],... %
               'parent',          [],... %Matrix containing each node's parent
               'links',           [],... %Matrix that holds the nodes' connections to each other
               'costs',           [],... %Matrix containing each connection's cost
               'probabilities',   [] ... %Matrix containing the probability for each connection
               ); 
 end

