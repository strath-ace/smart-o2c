function [checktot] = MyCreatedNodeCheck(Inputs, Attributes, ListNodes, parent)
% This function does the further checking on whether the node is valid
%
% Inputs:
% * Inputs    : The initialized inputs
% * newNode   : The structure containing the current node
% * ListNodes : The structure containing all the nodes
%
% Outputs: 
% * checktot  : The confirmation whether the node is a valid child. This
%               should be 1 if it is, and 0 otherwise

% Author: Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk

%Retreive the boundaries set for the node validity checks
CheckBounds = Inputs.NodeCheckBoundaries;

%Max dV_dep check
if (Attributes.tof==Attributes.tof_tot)
    check1 = (norm(Attributes.dV_dep) <= CheckBounds(1));
    %check2 = (Attributes.tof*86400*1e-7<= CheckBounds(1));
else
    check1 = (norm(Attributes.dV_dep) <= CheckBounds(2));
    %check2 = (Attributes.tof*86400*1e-7<= CheckBounds(2));
end

%Minimum perihelion check
check3 = ((Attributes.kep_trans.a*(1-Attributes.kep_trans.e)) >= CheckBounds(3));

%Low-Thrust Check
check4 = (Attributes.tof*86400*1e-7 >= CheckBounds(4)*(norm(Attributes.dV_dep)));

%Waiting time check
check5 = Attributes.t_dep-ListNodes.(parent).attributes.t_arr < CheckBounds(5);

%dV so far check
check6 = Attributes.dV_tot < CheckBounds(6);

%Edelbaum check
deltaV_edelbaum = sqrt( norm3(Attributes.v_dep)^2 + norm3(Attributes.lambertV_final)^2 - 2*norm3(Attributes.v_dep)*norm3(Attributes.lambertV_final));
check7 = Attributes.tof*86400*1e-7>=deltaV_edelbaum;

checktot = check1*check3*check4*check5*check6*check7;

end


