function [checktot] = MyCreatedNodeCheckMainBelt(Inputs, Attributes, ListNodes, parent)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%Retreive the boundaries set for the node validity checks
CheckBounds = Inputs.NodeCheckBoundaries;

%1 = check ok 
%0 = check fail (boundary exceeded)

%Max dV_dep check
check1 = (norm(Attributes.dV_dep) <= CheckBounds(1));

%Minimum perihelion check
check2 = ((Attributes.kep_trans.a*(1-Attributes.kep_trans.e)) >= CheckBounds(2));

%Low-Thrust Check
check3 = ((Attributes.t_arr-ListNodes.(parent).attributes.t_arr -20)*86400*1e-7 >= CheckBounds(3)*(norm(Attributes.dV_dep)));

%Low-Thrust Max dV check
%check4 = (Attributes.tof*86400*1e-7<= CheckBounds(5));

%Waiting time check
check5 = Attributes.t_dep-ListNodes.(parent).attributes.t_arr < CheckBounds(4);

%dV so far check
check6 = Attributes.dV_tot < CheckBounds(5);

checktot = check1*check2*check3*check5*check6;


end

