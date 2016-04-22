function [Attributes] = SetNodeAttributes(Inputs,characteristics)
% This function creates an object using the function specified in the
% options and sets attributes that identify the node using the
% characteristics
%s
% Inputs:
% * Inputs            : Structure containing the PhysarumSolver inputs
% * characteristics    : Vector containing the attributes to be assigned
%
% Outputs: 
% * Attributes   : Object containing the node's attributes
%
% Author: Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk
    

Attributes = Inputs.NodeAttributes();
attributenames = fieldnames(Attributes);
for i = 1:length(characteristics)
    attributeindex = Inputs.AttributeIDIndex(i);
    Attributes.(char(attributenames(attributeindex))) = characteristics(i);
end

end

