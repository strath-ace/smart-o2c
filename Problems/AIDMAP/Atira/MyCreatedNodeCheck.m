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
if (Attributes.tof==Attributes.tof_tot)
    check1 = (norm(Attributes.dV_dep) <= CheckBounds(1));
    %check2 = (Attributes.tof*86400*1e-7<= CheckBounds(1));
else
    check1 = (norm(Attributes.dV_dep) <= CheckBounds(2));
    %check2 = (Attributes.tof*86400*1e-7<= CheckBounds(2));
end

%Minimum perihelion check
check3 = ((Attributes.kep_trans.a*(1-Attributes.kep_trans.e)) >= CheckBounds(3));

%Low-Thrust Check (Edelbaum)
check4 = (Attributes.tof*86400*1e-7 >= CheckBounds(4)*(norm(Attributes.dV_dep)));

%Waiting time check
check5 = Attributes.t_dep-ListNodes.(parent).attributes.t_arr < CheckBounds(5);

%dV so far check
check6 = Attributes.dV_tot < CheckBounds(6);

checktot = check1*check3*check4*check5*check6;


% 
%  % Edelbaum's feasibility Check 
%             dV_edelbaum = sqrt( norm3(Attributes.dV_dep)^2 + norm3(Attributes.lambertV_final)^2 - 2*norm3(Attributes.v_dep)*norm3(Attributes.lambertV_final)); % [km/S]
%             deltaV_max      = ( global_config.low_trust_dVmax_coef * deltaV_Departure); % [km/S]
%             
%             if deltaV_edelbaum > deltaV_max
%                 deltaV_limit   = deltaV_edelbaum;
%             else 
%                 deltaV_limit   = deltaV_max;
%             end
%             
%             % Compute maximum dv provided by the Low-Thrust engine [m/s]
%             dvlowthrust = global_config.low_thrust_acc * curr_tof * 86400;
%         
%             % Sanity Check: Low-Thrust Constrain            
%             if dvlowthrust < ( deltaV_limit * 1000) % [m/s]                      
%                 continue;
%             end         
% 

end


