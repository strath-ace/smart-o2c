function [newNode] = MyCreatedNodeCheck(Inputs, newNode, ListNodes)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


%1 = check ok 
%0 = check fail (boundary exceeded)

%Max dV_dep check
check1 = (norm(newNode.attributes.dV_dep) <= 3);

%Minimum perihelion check
check2 = ((newNode.attributes.kep_trans.a*(1-newNode.attributes.kep_trans.e)) >= 0.31);

%Low-Thrust Check
check3 = (newNode.attributes.tof*86400*1e-7 >= 2*(norm(newNode.attributes.dV_dep)));

%Waiting time check
check4 = newNode.attributes.t_dep-ListNodes.(newNode.parent).attributes.t_arr < 2*365;

% prevNode = ListNodes.(newNode.parent);
% prev_orbit = prevNode.attributes.kep_trans;
% 
% [~, departure_v] = StardustTool.CartesianElementsAt(prev_orbit,newNode.attributes.t_dep);
% p = newNode.attributes.kep_trans.a*(1-newNode.attributes.kep_trans.e^2);
% V0 = norm(departure_v);
% Vf = norm(newNode.attributes.lambertV_final);
% 
% dVEdelbaum = sqrt(V0^2-2*V0*Vf*cos(pi/2*newNode.attributes.kep_trans.i*pi/180)+Vf^2);
% 
% check5 = (newNode.attributes.tof*86400*1e-7 >= 2*dVEdelbaum);

checktot = check1*check2*check3*check4;
%checktot = check1*check2*check3*check4*check5;




if checktot ==1
    newNode.length = newNode.length;
else
    newNode.length = inf;
end

end

