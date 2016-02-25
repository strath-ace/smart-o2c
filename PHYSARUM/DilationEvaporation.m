function [Nodes] = DilationEvaporation(Nodes,Inputs,WalkedPaths)
% This function updates the radius of links walked by the agents with the
% dilation and evaporation mechanics.
%
% Inputs:
% * Nodes       : Structure containing the graph
% * Inputs      : Structure containing the PhysarumSolver inputs
% * WalkedPaths : The paths walked by the agent (matrix of the same size as the Nodes.links matrix,
%                 with a 1 for each path walked 
%
% Outputs: 
% * Nodes       : Nodes structure where the radii have been updated with
%                 the dilation and evaporation
%
% Author: Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk

%Calculate L_tot (the cost)
L_tot = sum(sum(WalkedPaths.*Nodes.lengths));

%Update radii with dilation
Nodes.radius = Nodes.radius + Inputs.LinearDilationCoefficient.*WalkedPaths.*Nodes.radius/L_tot;

%Update radii with evaporation
Nodes.radius = Nodes.radius - Inputs.EvaporationCoefficient.*WalkedPaths.*Nodes.radius;

%Prevent radii from becoming too large (exploding) or small (closing)
Nodes.radius(Nodes.radius > Inputs.MaximumRadius) = Inputs.MaximumRadius;
Nodes.radius(Nodes.radius < Inputs.MinimumRadius) = Inputs.MinimumRadius;

end

