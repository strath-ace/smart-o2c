function [newNode] = MyCreatedNodeCheck(Inputs, newNode, ListNodes)
% This function does the further checking on whether the node is valid
%
% Inputs:
% * Inputs   : The initialized inputs
% * newNode   : The structure containing the current node
% * ListNodes : The structure containing all the nodes
%
% Outputs: 
% * newNode   : The updated node structure

% Author: Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk


%Retreive the boundaries set for the node validity checks
CheckBounds = Inputs.NodeCheckBoundaries;

%1 = check ok 
%0 = check fail (boundary exceeded)

%Max dV_dep check
check1 = (norm(newNode.attributes.dV_dep) <= CheckBounds(1));

%Minimum perihelion check
check2 = ((newNode.attributes.kep_trans.a*(1-newNode.attributes.kep_trans.e)) >= CheckBounds(2));

%Low-Thrust Check
check3 = (newNode.attributes.tof*86400*1e-7 >= CheckBounds(3)*(norm(newNode.attributes.dV_dep)));

%Waiting time check
check4 = newNode.attributes.t_dep-ListNodes.(newNode.parent).attributes.t_arr < CheckBounds(4);

checktot = check1*check2*check3*check4;

if checktot ==1
    newNode.length = newNode.length;
else
    newNode.length = inf;
end

end

