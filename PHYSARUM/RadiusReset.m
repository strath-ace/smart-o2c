function [ output_args ] = RadiusReset(Inputs,ListNodes)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

nodenames = char(fieldnames(ListNodes));
for i = 1:length(nodenames)
    ListNodes.(nodenames).radius = Inputs.StartingRadius*ones(1, length(ListNodes.(nodenames).radius));
    for j = 1:length(ListNodes.(nodenames).radius)
        

end

