function [newnode_ID] = ChooseNode(possnodes)
% This function is used to choose a new node for ramification
%
% Inputs:
% * possnodes   : Vector containing the node_IDs of all the possible nodes
% 
% Outputs:
% * newnode_ID        : The node_ID chosen
%
% Author: Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk


%Choose a random decision (node_ID)
randdecision = randi([1 length(possnodes)]);
newnode_ID = char(possnodes(randdecision));


end

