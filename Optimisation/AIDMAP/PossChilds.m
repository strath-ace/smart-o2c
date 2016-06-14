function [posschildren] = PossChilds(targets,sets)
% This script creates a cell array containing all the unique IDs (UIDs) the
% Physarum algorithm can choose from. 
%
% Inputs:
% * targets    : A cell array with the names of the targets
% * sets       : A structure with the possible values each attribute
%                shown in the unique ID can take.
%
% Outputs: 
% * posschildren  : A cell array with all the UIDs the algorithm can
%                choose from
%
% Author: Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk

%Retrieve the possible attribute values & create the posschildren vector
attribvalues = struct2cell(sets);
posschildren =[];

%Loop over all the targets
for i = 1:length(targets)
    
    %Set for convenience the current target name as a variable
    targetname = targets(i);
    
    %Set the first column as the target's index
    possibleattributes{i,1} = i;
    
    %Loop over all the attributes
    for j = 1:length(attribvalues)
       
        %Add the attribute values to the posssibleattributes array       
        possibleattributes{i,j+1} = 1:1:length(cell2mat(attribvalues{j}(i)));
    end
           
    %Determine all possible combinations of the possible characteristics
    %for this target
    possiblecombinations = combvec(possibleattributes{i,:})';
    
    
    %%%Create UID%%%
    %As only underscores can be used in field names, the UID is made up as follows:
    %3 underscores seperate the parent's section of the ID and the child's
    %2 underscores define the difference between the target and the attributes
    %1 underscores define the difference between two attributes

    %First, combine all attributes into 1 string, where the space between
    %the attributes is filled with underscores. The goal is to separate
    %them by exactly 1 underscore, but this is a necessary step to do so
    for j = 1:length(possiblecombinations)
        temp = possiblecombinations(j,:);
        possattribstr{j} = strrep(num2str(temp(:)'),' ','_');
    end    
    
    %As MATLAB sometimes sometimes replaces multiple spaces by underscores,
    %a while loop is needed to correct the strings such that they only
    %contain single consecutive underscores.
    while ~isempty(cell2mat(strfind(possattribstr,'__')))
        possattribstr = strrep(possattribstr,'__','_');
    end
        
    %Combine the target name and the attributes into the UID and add them
    %to a vector that will contain all the possible children
    posschildren = [posschildren strcat(targetname,'__',possattribstr)];
    clear possattribstr
end

end

