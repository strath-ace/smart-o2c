function [checktot] = MyCreatedNodeCheck(Inputs, Attributes, ListNodes, parent)
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%-----------Copyright (C) 2018 University of Strathclyde and Authors-----------
%
%
%
%% MyCreatedNodeCheck: This function does the further checking on whether the node is valid
% 
%% Inputs:
% * Inputs      : The initialised inputs structure
% * Attributes  : The structure containing the problem-specific attributes of the node to be checked
% * ListNodes   : The structure containing all the nodes
% * parent      : The parent of the node to be checked [string]
% 
%% Outputs: 
% * checktot  : The confirmation whether the node is a valid child. This
%               should be 1 if it is, and 0 otherwise 
% 
%% Author(s): Aram Vroom (2016)
% Email:  aram.vroom@strath.ac.uk

% Retreive the boundaries set for the node validity checks
CheckBounds = Inputs.NodeCheckBoundaries;

% Max dV_dep check
if (Attributes.tof==Attributes.tof_tot)
    check1 = (norm(Attributes.dV_dep) <= CheckBounds(1));
else
    check1 = (norm(Attributes.dV_dep) <= CheckBounds(2));
end

% Minimum perihelion check
check2 = ((Attributes.kep_trans.a*(1-Attributes.kep_trans.e)) >= CheckBounds(3));

% Low-Thrust Check
check3 = (Attributes.tof*86400*1e-7 >= CheckBounds(4)*(norm(Attributes.dV_dep)));

% Waiting time check
check4 = Attributes.t_dep-ListNodes.(parent).attributes.t_arr < CheckBounds(5);

% dV so far check
check5 = Attributes.dV_tot < CheckBounds(6);

% Edelbaum check
deltaV_edelbaum = sqrt( norm3(Attributes.v_dep)^2 + norm3(Attributes.lambertV_final)^2 - 2*norm3(Attributes.v_dep)*norm3(Attributes.lambertV_final));
check6 = Attributes.tof*86400*1e-7>=deltaV_edelbaum;

% Combine all performed checks
checktot = check1*check3*check4*check5*check6;

end


