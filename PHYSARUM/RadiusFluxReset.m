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
    %Set the vein radii to their default value
    ListNodes.(char(nodenames(i))).radius = Inputs.StartingRadius*ones(1, length(ListNodes.(char(nodenames(i))).radius));
    
    %Loop over all the children
    for j = 1:length(ListNodes.(char(nodenames(i))).children)
        
        %Recalculate the fluxes using the previously updated radius
        ListNodes.(char(nodenames(i))).fluxes(j) = CalculateFlux(Inputs, ListNodes, ListNodes.(char(nodenames(i))).children(j));
    end
end
        

end

