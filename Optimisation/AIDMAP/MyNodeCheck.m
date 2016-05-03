function [validflag] = MyNodeCheck(ListNodes,newnode_ID,currentNode,generatednodes)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%Check 2
check2 = isempty(strmatch(newnode_ID, fields(generatednodes), 'exact'));

%Check 3 - ToF check, assuming ToF = 1st attribute in ID & Ta = 2nd
%attribute

%Extract child node
temp = strsplit(newnode_ID,'____');
childnode = temp{2};

temp = strsplit(childnode,'___');
attribs = temp{2};
temp = strsplit(char(attribs),'__');
chosentof = str2double(strrep(temp(1),'_','.'));
chosent_arr = str2double(strrep(temp(2),'_','.'));

parentt_arr = ListNodes.(char(currentNode)).attributes.t_arr;

check3 = (chosent_arr - chosentof > parentt_arr);




validflag = check2*check3;            

end

