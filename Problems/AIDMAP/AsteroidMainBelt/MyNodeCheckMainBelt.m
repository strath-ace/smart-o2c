function [validflag] = MyNodeCheckMainBelt(Inputs,ListNodes,newnode_ID,currentNode,generatednodes)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

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

