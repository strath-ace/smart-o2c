function [newNode] = MyCreatedNodeCheck(Inputs, newNode, ListNodes)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


%1 = check ok 
%0 = check fail (boundary exceeded)
check1 = (norm(newNode.attributes.dV_dep) <= 3);
check2 = ((newNode.attributes.kep_trans.a*(1-newNode.attributes.kep_trans.e)) >= 0.31);
check3 = (newNode.attributes.tof*86400*1e-7 >= 2*(norm(newNode.attributes.dV_dep)));


prevNode = ListNodes.(newNode.parent);
prev_orbit = prevNode.attributes.kep_trans;

[~, departure_v] = StardustTool.CartesianElementsAt(prev_orbit,newNode.attributes.t_dep);
p = newNode.attributes.kep_trans.a*(1-newNode.attributes.kep_trans.e^2);
V0 = norm(departure_v);
Vf = norm(newNode.attributes.lambertV_final);

dVEdelbaum = sqrt(V0^2-2*V0*Vf*cos(pi/2*newNode.attributes.kep_trans.i*pi/180)+Vf^2);

check4 = (newNode.attributes.tof*86400*1e-7 >= dVEdelbaum);
checktot = check1*check2*check3*check4;

if checktot ==1
    newNode.length = newNode.length;
else
    newNode.length = inf;
end

end

