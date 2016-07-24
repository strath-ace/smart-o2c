function [Flux] = CalculateFlux(Inputs, Node)
%% CalculateFlux: This function calculates the flux for a specified node
% 
%% Inputs:
% * Inputs      : Structure containing the PhysarumSolver inputs
% * Node        : Structure that cotnains the node and its attributes
% 
%% Outputs:
% * Flux        : The flux for the node [real number]
% 
%% Author: Aram Vroom (2016)
% Email:  aram.vroom@strath.ac.uk

% Retreive the link's radius
radius = Node.radius;

% Find the link's length
length = Node.length;

% Calculate the flux
Flux = pi*radius^4/(8*Inputs.Viscosity)*1/length;
end

