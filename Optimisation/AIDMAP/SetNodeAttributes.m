function [Attributes] = SetNodeAttributes(Inputs, Parent, cityname, attributes)
%% SetNodeAttributes: This function creates an object using the function specified in the options and sets attributes that identify the node using the characteristics
%
%% Inputs:
% * Inputs            : Structure containing the PhysarumSolver inputs
% * Parent            : String containing the name of the node's parent
% * cityname          : Name of the current city
% * attributes        : Vector containing the attributes to be assigned
%
%% Outputs: 
% * Attributes   : Object containing the node's attributes
%
%% Author: Aram Vroom (2016)
% Email:  aram.vroom@strath.ac.uk
    
%Retrieve the node attributes defined by the user
Attributes = Inputs.MyNodeAttributes();
SetNames = fieldnames(Inputs.Sets);

%Obtain their names
attributenames = fieldnames(Attributes);

%Check if the current node is the root
if strcmp(cityname,Inputs.RootName)
    %Loop over the attributes
    for i = 1:length(attributes)

        %Find the index of the attributes set by the user
        attributeindex = Inputs.AttributeIDIndex(i);

        %Set the value of these attributes as defined
        Attributes.(char(attributenames(attributeindex))) = Inputs.RootAttrib(i);
    end
else
    %If the current node is not the root, start from i = 2, as i = 1 is the
    %city's index (which is not the case for the root)
    for i = 2:length(attributes)
        
        %Find the index of the attributes set by the user
        attributeindex = Inputs.AttributeIDIndex(i-1);
        
        %Set the value of these attributes as defined
        Attributes.(char(attributenames(attributeindex))) = Inputs.Sets.(char(SetNames(i-1))){attributes(1)}(attributes(i));
        
    end
end
        
        
%Calculate the additional attibutes
Attributes = Inputs.MyAttributeCalcFile(Inputs, Parent, cityname, Attributes);

end

