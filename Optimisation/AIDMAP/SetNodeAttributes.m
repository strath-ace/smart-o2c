function [Attributes] = SetNodeAttributes(Inputs, Parent,targetname,attributes)
% This function creates an object using the function specified in the
% options and sets attributes that identify the node using the
% characteristics
%s
% Inputs:
% * Inputs            : Structure containing the PhysarumSolver inputs
% * attributes        : Vector containing the attributes to be assigned
%
% Outputs: 
% * Attributes   : Object containing the node's attributes
%
% Author: Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk
    
%Retrieve the node attributes defined by the user
Attributes = Inputs.NodeAttributes();

%Obtain their names
attributenames = fieldnames(Attributes);

%Loop over the attributes
for i = 1:length(attributes)
    
    %Find the index of the attributes set by the user
    attributeindex = Inputs.AttributeIDIndex(i);
    
    %Set the value of these attributes as defined
    Attributes.(char(attributenames(attributeindex))) = attributes(i);
end

Attributes = Inputs.MyAttributeCalcFile(Inputs, Parent, targetname, Attributes);

end

