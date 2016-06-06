function [validflag] = MyNodeCheck(Inputs,ListNodes,newnode_ID,currentNode,generatednodes)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%Check 2
check2 = isempty(strmatch(newnode_ID, fields(generatednodes), 'exact'));

%Check 3 - ToF check, assuming ToF = 1st attribute in ID & Ta = 2nd
%attribute

%Extract child node
temp = strsplit(newnode_ID,'___');
childnode = temp{2};

temp = strsplit(childnode,'__');
attribs = temp{2};
temp = strsplit(char(attribs),'_');

asteroidindex = sscanf(temp{1},'%i');
tofindex = sscanf(temp{2},'%i');
t_arrindex = sscanf(temp{3},'%i');

chosentof = Inputs.Sets.tof{asteroidindex}(tofindex);
chosent_arr = Inputs.Sets.epochsnode{asteroidindex}(t_arrindex);

parentt_arr = ListNodes.(char(currentNode)).attributes.t_arr;


check3 = (chosent_arr - chosentof > parentt_arr);




validflag = check2*check3;            

end

