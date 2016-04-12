function [Cost] = CostFunction(fromNode, toNode)
% This function calculates the cost of a certain connection.
% It can be altered such that it is applicable to the problem at hand.
%
% Inputs:
% * fromNode  : The node from which the cost is calculated
% * Inputs    : The node to which the cost should be calculated
%
% Outputs: 
% * Cost   : The cost of the link between the two nodes
%
% Author: Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk
    

%Calculate cost to change orbital elements with orbit characterstics such
%as ToF set. Currently simple formula to test functionality.
Cost = (toNode.characteristics(1)-fromNode.characteristics(1))^2;