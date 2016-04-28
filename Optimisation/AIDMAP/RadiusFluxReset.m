function [ListNodes] = RadiusFluxReset(Inputs, ListNodes)
% This function handles the restarting of the Physarum algorithm by
% resetting the radii to their default values and updating the fluxes accordingly.
%
% Inputs:
% * Inputs   : Structure containing the PhysarumSolver inputs
% * ListNodes       : Structure containing the graph
%
% Outputs: 
% * ListNodes       : Structure containing the updated graph
%
% Author: Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk

%Retrieve the node IDs
nodenames = fieldnames(ListNodes);

%Loop over all the nodes
for i = 1:length(nodenames)
    
    if ~isempty(ListNodes.(char(nodenames(i))).parent)  %Ignore root
    
    %Set the vein radii to their default value
    ListNodes.(char(nodenames(i))).radius = Inputs.StartingRadius*ones(1, length(ListNodes.(char(nodenames(i))).radius));
    
    %Recalculate the fluxes using the previously updated radius
    ListNodes.(char(nodenames(i))).flux = CalculateFlux(Inputs, ListNodes.(char(nodenames(i))));
    
    end

end
        

end

