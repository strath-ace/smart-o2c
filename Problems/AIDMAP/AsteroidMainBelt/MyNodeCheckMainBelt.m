function [validflag] = MyNodeCheckMainBelt(Inputs, ListNodes, newnode_ID, currentNode)
%% MyNodeCheckMainBelt: This function checks the validity of the node to be created using solely the unique ID 
% 
%% Inputs:
% * Inputs      : The initialised inputs structure
% * ListNodes   : The structure containing all the nodes
% * newNode     : The node ID of the node to be created [string]
% * currentNode : The ID of the current node [string]
% 
%% Outputs: 
% * validflag   : The flag that indicates whether the node to be created is
%                 valid. This should be 1 if the node is valid and 0
%                 otherwise
% 
%% Author: Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk

% Check 1 - ToF check: confirm that the required departure date is not
% before the current time

% Obtain attribute indices
temp = strsplit(newnode_ID, '_');

asteroidindex = sscanf(temp{end-2}, '%f');
tofindex = sscanf(temp{end-1}, '%f');
t_arrindex = sscanf(temp{end}, '%f');

% Obtain the time of flight and arrival time of the node
chosentof = Inputs.Sets.tof{asteroidindex}(tofindex);
chosent_arr = Inputs.Sets.epochsnode{asteroidindex}(t_arrindex);

% Retrieve the arrival date of the parent
parentt_arr = ListNodes.(char(currentNode)).attributes.t_arr;

% Check the arrival time - the time of flight is larger than the arrival
% time of the parent
check1 = (chosent_arr - chosentof > parentt_arr);

% Set the flag
validflag = check1;                 
                  

end

