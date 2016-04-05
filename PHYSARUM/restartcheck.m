function [restartflag] = restartcheck(Nodes,currentnode)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

restartflag = 0;

if isempty(Nodes.ListNodes.(char(currentnode)).possibledecisions)
    %disp(strcat(currentnode, ' Has no More Targets'))
    restartflag = 1;
    return
end


end

