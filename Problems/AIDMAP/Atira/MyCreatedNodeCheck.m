function [checktot] = MyCreatedNodeCheck(Inputs, Attributes, ListNodes, parent)
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
check1 = (norm(Attributes.dV_dep) <= CheckBounds(1));

%Minimum perihelion check
check2 = ((Attributes.kep_trans.a*(1-Attributes.kep_trans.e)) >= CheckBounds(2));

%Low-Thrust Check
check3 = (Attributes.tof*86400*1e-7 >= CheckBounds(3)*(norm(Attributes.dV_dep)));

%Waiting time check
check4 = Attributes.t_dep-ListNodes.(parent).attributes.t_arr < CheckBounds(4);

checktot = check1*check2*check3*check4;


end

