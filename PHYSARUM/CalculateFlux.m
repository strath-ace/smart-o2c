function [Flux] = CalculateFlux(Inputs,ListNodes,node)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

node = char(node);

parent = ListNodes.(node).parent;
radius = ListNodes.(node).radius
length = ListNodes.(node).length;

Flux = pi*radius^4/(8*Inputs.Viscosity)*1/length;
end

