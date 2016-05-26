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

%Loop over all the targets

attribvalues = struct2cell(sets);

for i = 1:length(targets)
    
    %Set for convenience the current target name as a variable
    targetname = targets(i);
    
    %Loop over all the attributes
    for j = 1:length(attribvalues)
        
        %Determine whether the currently evaluated attribute has a
        %matrix or vector as input
        vectorflag = isvector(attribvalues{j});
        onevalueflag = (length(attribvalues{j}(1,:))==1);
        cellflag = iscell(attribvalues{j});
        
        possibleattributes{i,j} = cell2mat(attribvalues{j}(i));
    end
           
    %Determine all possible combinations of the possible characteristics
    %for this target
    possiblecombinations = combvec(possibleattributes{i,:})';
    
    
    %%%Create UID%%%
    %As only underscores can be used in field names, the UID is made up as follows:
    %3 underscores define the difference between the target and the
    %  attributes
    %2 underscores define the difference between two attributes
    %1 underscore denotes the location of the decimal point within an
    %  attribute
    
    %First, combine all attributes into 1 string, where the different
    %attributes are initially seperated by 1 underscore. The goal is
    %to have each attribute separated by two, but this is a neccessary
    %the initial step to do so.
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
    
    %Once all the attributes are separated by exactly 1 underscore, 
    %replace them all by 2 underscores to reach the final result.
    possattribstr = strrep(possattribstr,'_','__');
    
    %Combine the target name and the attributes into the UID by
    %looping over this target's possible characteristics and concatenating
    %the two strings.
    for j = 1:length(possattribstr)
        possdecattribvec(i,j) = strcat(targetname,'___',possattribstr{j});
        
        %As dots are not possible in field names, use 1 underscore to
        %denote the location of the decimal point
        possdecattribvec(i,j) = strrep(possdecattribvec(i,j),'.','_');
    end
    clear possattribstr
end


%Reshape the possdecattribvec variable such that it becomes a cell array 
%with a single row.
posschildrenvec = reshape(possdecattribvec, [1, numel(possdecattribvec)]);

%Remove empty cells that are the result of more possible characteristics
%for certain targets
posschildren = posschildrenvec(~cellfun('isempty',posschildrenvec));
end

