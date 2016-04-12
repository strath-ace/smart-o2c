function [ListNodes] = RadiusFluxReset(Inputs,ListNodes)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

nodenames = fieldnames(ListNodes);
for i = 1:length(nodenames)
    ListNodes.(char(nodenames(i))).radius = Inputs.StartingRadius*ones(1, length(ListNodes.(char(nodenames(i))).radius));
    for j = 1:length(ListNodes.(char(nodenames(i))).radius)
        ListNodes.(char(nodenames(i))).fluxes(j) = CalculateFlux(Inputs, ListNodes, ListNodes.(char(nodenames(i))).children(j));
    end
end
        

end

