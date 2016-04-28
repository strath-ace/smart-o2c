function [posschildren] = PossChilds(targets,sets)
% This script creates a cell array containing all the unique IDs (UIDs) the
% Physarum algorithm can choose from. 
%
% Inputs:
% * targets    : A cell array with the names of the targets
% * charvalues : A cell array with the possible values each characteristic
%                shown in the unique ID can take. If the possible values for 
%                a characteristic are the same for each target, a vector can 
%                be used. Otherwise, a matrix should be the input. 
%
% Outputs: 
% * posschildren  : A cell array with all the UIDs the algorithm can
%                choose from
%
% Author: Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk

%Loop over all the targets

charvalues = struct2cell(sets);

for i = 1:length(targets)
    
    %Set for convenience the current target name as a variable
    targetname = targets(i);
    
    %Loop over all the characteristics
    for j = 1:length(charvalues)
        
        %Determine whether the currently evaluated characterstic has a
        %matrix or vector as input
        vectorflag = isvector(charvalues{j});
        onevalueflag = (length(charvalues{j}(1,:))==1);
        cellflag = iscell(charvalues{j});
        
        possiblecharacteristics{i,j} = cell2mat(charvalues{j}(i));
    end
           
    %Determine all possible combinations of the possible characteristics
    %for this target
    possiblecombinations = combvec(possiblecharacteristics{i,:})';
    
    
    %%%Create UID%%%
    %As only underscores can be used in field names, the UID is made up as follows:
    %3 underscores define the difference between the target and the
    %  characteristics
    %2 underscores define the difference between two characteristics
    %1 underscore denotes the location of the decimal point within a
    %  characteristic
    
    %First, combine all characteristics into 1 string, where the different
    %characteristics are initially seperated by 1 underscore. The goal is
    %to have each characteristic separated by two, but this is a neccessary
    %the initial step to do so.
    for j = 1:length(possiblecombinations)
        temp = possiblecombinations(j,:);
        posscharrstr{j} = strrep(num2str(temp(:)'),' ','_');
    end
    
    
    %As MATLAB sometimes sometimes replaces multiple spaces by underscores,
    %a while loop is needed to correct the strings such that they only
    %contain single consecutive underscores.
    while ~isempty(cell2mat(strfind(posscharrstr,'__')))
        posscharrstr = strrep(posscharrstr,'__','_');
    end
    
    %Once all the characteristics are separated by exactly 1 underscore, 
    %replace them all by 2 underscores to reach the final result.
    posscharrstr = strrep(posscharrstr,'_','__');
    
    %Combine the target name and the characteristics into the UID by
    %looping over this target's possible characteristics and concatenating
    %the two strings.
    for j = 1:length(posscharrstr)
        possdeccharvec(i,j) = strcat(targetname,'___',posscharrstr{j});
        
        %As dots are not possible in field names, use 1 underscore to
        %denote the location of the decimal point
        possdeccharvec(i,j) = strrep(possdeccharvec(i,j),'.','_');
    end
end


%Reshape the possdeccharvec variable such that it becomes a cell array 
%with a single row.
posschildrenvec = reshape(possdeccharvec, [1, numel(possdeccharvec)]);

%Remove empty cells that are the result of more possible characteristics
%for certain targets
posschildren = posschildrenvec(~cellfun('isempty',posschildrenvec));
end

