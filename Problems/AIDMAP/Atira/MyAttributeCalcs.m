function [Attributes] = MyAttributeCalcs(Inputs, Parent, targetname, Attributes)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%Check if the target is in the Asteroid class
if ismember(targetname,fieldnames(Asteroids))
    
    %Retrieve the cartesian coordiantes of the target at the arrival epoch
    Attributes.r_arr = StardustTool.CartesianElementsAt(Asteroids.(targetname),Attributes.t_arr);
    
    %Find the departure time by subtracting the ToF from the arrival time
    Attributes.t_dep = Attributes.t_arr - Attributes.tof;
    
    %Check if the node has a parent (excludes the root)
    if ~isempty(Parent)
        
        %Add the ToF for this transfer to the total ToF of the parent
        %to find the new total ToF
        Attributes.tof_tot = Parent.attributes.tof_tot + Attributes.tof;
    else
        
        %If no parent has been set, the total ToF is the current ToF
        Attributes.tof_tot = Attributes.tof;

        %The current keplerian orbit is then assumed to be that of the root
        kep = Asteroids.(targetname).getKeplerianElements;
        Attributes.kep_trans = [kep.a kep.e kep.i kep.OM kep.W kep.M0 kep.t0];
    end 
    
else
    
    %If the target is not in the Asteroids list, set the r_arr, tof_tot and
    %t_dep to default values.
    Attributes.r_arr = zeros(1,3);
    Attributes.tof_tot = 0;
    Attributes.t_dep = 0;
end


end
