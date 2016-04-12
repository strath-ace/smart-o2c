function [Flux] = CalculateFlux(Inputs,ListNodes,childnode)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

childnode = char(childnode);

parent = ListNodes.(childnode).parent;
childindex = find(strcmp(ListNodes.(parent).children,childnode)==1);
radius = ListNodes.(parent).radius(childindex);
length = ListNodes.(parent).lengths(childindex);

Flux = pi*radius^4/(8*Inputs.Viscosity)*1/length;
end

