function [Nodes] = CreateGraph()
% This function creates the structure describing the graph
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

Nodes = struct('identifiers',        [],... % Identifier of each node
               'parents',            [],... % Matrix containing each node's parent
               'links',              [],... % Matrix that holds the nodes' connections to each other
               'radii',              [],... % The radius of each connection
               'pressure_gradients', [],... % The pressure gradient over each connection
               'lengths',            [],... % The length of each connection
               'fluxes',             [],... % Matrix containing each connection's flux
               'probabilities',      [] ... % Matrix containing the probability for each connection
               ); 
 end

