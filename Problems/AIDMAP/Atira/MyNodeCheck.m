function [validflag] = MyNodeCheck(ListNodes,newnode_ID,currentNode,generatednodes)
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
check2 = isempty(strmatch(newnode_ID, fields(generatednodes), 'exact'));

%Check 3 - ToF check, assuming ToF = 1st attribute in ID & Ta = 2nd
%attribute

%Extract child node
temp = strsplit(newnode_ID,'____');
childnode = temp{2};

%Find the attributes put in the ID
temp = strsplit(childnode,'___');
attribs = temp{2};

%Find the chosen ToF & t_arr
temp = strsplit(char(attribs),'__');
chosentof = str2double(strrep(temp(1),'_','.'));
chosent_arr = str2double(strrep(temp(2),'_','.'));

%Obtain the parent's t_arr
parentt_arr = ListNodes.(char(currentNode)).attributes.t_arr;

%Check departure time isn't before the current time
check3 = (chosent_arr - chosentof > parentt_arr);


validflag = check2*check3;            

end

