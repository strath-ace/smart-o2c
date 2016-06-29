function [validflag] = MyNodeCheck(Inputs, ListNodes,newnode_ID,currentNode,generatednodes)
% This function checks the validity of the node to be created using the unique ID 
%
% Inputs:
% * Inputs      : The initialized inputs
% * ListNodes   : The structure containing all the nodes
% * newNode     : The node ID of the node to be created
% * currentNode : the ID of the current node
%
% Outputs: 
% * validflag   : The flag that indicates whether the node to be created is
%                 valid

% Author: Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk

%Check 2
check2 = isempty(find(strcmp(newnode_ID, fields(generatednodes)), 1));

%Check 3 - ToF check, assuming ToF = 1st attribute in ID & Ta = 2nd
%attribute

%Obtain attribute indices
temp = strsplit(newnode_ID,'_');

asteroidindex = sscanf(temp{end-2}, '%f');
tofindex = sscanf(temp{end-1}, '%f');
t_arrindex = sscanf(temp{end}, '%f');

chosentof = Inputs.Sets.tof{asteroidindex}(tofindex);
chosent_arr = Inputs.Sets.epochsnode{asteroidindex}(t_arrindex);

parentt_arr = ListNodes.(char(currentNode)).attributes.t_arr;


check3 = (chosent_arr - chosentof > parentt_arr);




validflag = check2*check3;                 
       

end

