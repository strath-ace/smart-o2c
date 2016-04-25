function [possnodes] = PossNodes(targets,charvalues)
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
% * possnodes  : A cell array with all the UIDs the algorithm can
%                choose from
%
% Author: Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk

%Loop over all the targets
for i = 1:length(targets)
    
    %Set for convenience the current target name as a variable
    targetname = targets(i);
    
    %Loop over all the characteristics
    for j = 1:length(charvalues)
        
        %Determine whether the currently evaluated characterstic has a
        %matrix or vector as input
        vectorflag = isvector(charvalues{j});
        cellflag = iscell(charvalues{j});
        
        if vectorflag
            
            %If the values of the characteristic have been defined as a
            %cell array, the value is different for each target and a cell
            %array was used due to the number of values of the
            %characteristic not being equal for each target (thus it not
            %being possible to create a matrix). As such, a different value
            %is taken for each target.
            if cellflag
                possiblecharacteristics{i,j} = cell2mat(charvalues{j}(i));
            else
                
            %Otherwise, the possible values for the characterstic are given as
            % a normal vector. The possible values are then the same for every target. 
            %As such, the same "charvalues" values are taken for every loop over i.
                possiblecharacteristics{i,j} = charvalues{j};
            end         
                 
        else
            
            %If the possible values are given as a matrix, these possible
            %values may differ per target, but the number of values the characteristic
            %is equal for each target. As such, the ith row is taken of
            %this charvalues cell.
            possiblecharacteristics{i,j} = charvalues{j}(i,:);
        end
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
possnodesvec = reshape(possdeccharvec, [1, numel(possdeccharvec)]);

%Remove empty cells that are the result of more possible characteristics
%for certain targets
possnodes = possnodesvec(~cellfun('isempty',possnodesvec));
end

