function [Cost] = MyCostFunction(fromNode, toNode)
% This function calculates the cost of a certain connection.
% It can be altered such that it is applicable to the problem at hand.
%
% Inputs:
% * fromNode  : The node from which the cost is calculated (structure)
% * toNode    : The node to which the cost is calculated (structure)
%
% Outputs: 
% * Cost   : The cost of the link between the two nodes
%
% Author: Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk

attributenames = fieldnames(toNode.characteristics);

%Calculate cost to change orbital elements with orbit characterstics such
%as ToF set. Currently simple formula to test functionality.
Cost = ((toNode.characteristics.t_arr - toNode.characteristics.tof) - fromNode.characteristics.t_arr)^5;
