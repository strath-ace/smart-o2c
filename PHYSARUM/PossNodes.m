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
                possiblecharacteristics{i,j} = cell2mat(charvalues{j}{i});
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
    
    %Combine each combination into 1 string, where the different
    %characteristics are seperated by an underscore
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
    
    %Combine the target name and the characteristics into the UID by
    %looping over this target's possible cahracteristics and concatenating
    %the two strings.
    for j = 1:length(posscharrstr)
        possdeccharvec(i,j) = strcat(targetname,'__',posscharrstr{j});
    end
end

%Reshape the possdeccharvec variable such that it becomes a cell array 
%with a single row.
possnodes = reshape(possdeccharvec, [1, numel(possdeccharvec)]);

end

