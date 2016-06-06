function [continueflag] = MyEndConditions(Inputs,Agents,agent)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

currentnode = Agents.(char(agent)).currentNode;

endconditions = Inputs.EndConditions;

%Check if final target reached
temp = strsplit(currentnode,'___');
temp = strsplit(temp{end},'__');
currenttarget = temp{1};

check1 = ~(sum(ismember(currenttarget, endconditions{1}))==length(currenttarget));

%Total dV
%check2 = ~(sum(Agents.(char(agent)).previouscosts)>=endconditions{2});

continueflag = check1;

%continueflag = check1*check2;
end

