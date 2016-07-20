function [continueflag] = MyEndConditions(Inputs,Agents,agent)
% This function confirms whether the final condtions have been reached
%
% Inputs:
% * Inputs        : Structure containing the PhysarumSolver inputs
% * Agents        : Structure containing the agents
% * agent         : The name of the agent currently being evaluated
%
% Outputs: 
% * continueflag  : The flag that confirms whether the algorithm has
%                   reached the final condition. This should be 1 if 
%                   it has, and 0 otherwise
%
% Author: Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk

%Retrieve the current node
currentnode = Agents.(char(agent)).currentNode;

%Obtain the final conditions
endconditions = Inputs.EndConditions;

%Check if final target reached
temp = strsplit(currentnode,'_');
currenttarget = temp{end-3};
check1 = ~(sum(ismember(currenttarget, endconditions{1}))==length(currenttarget));

%Determine the continueflag
continueflag = check1;

end
