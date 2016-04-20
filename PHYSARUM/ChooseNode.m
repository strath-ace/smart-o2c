function [newnode_ID] = ChooseNode(possnodes)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%Choose a random decision (node_ID)
randdecision = randi([1 length(possnodes)]);
newnode_ID = char(possnodes(randdecision));


end

