function [checktot] = MyCreatedNodeCheckMainBelt(Inputs, Attributes, ListNodes, parent)
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%-----------Copyright (C) 2016 University of Strathclyde-------------
%
%
%
%% MyCreatedNodeCheckMainBelt: This function does the further checking on whether the node is valid
% 
%% Inputs:
% * Inputs      : The initialised inputs structure
% * Attributes  : The structure containing the problem-specific attributes 
%                 of the node to be checked
% * ListNodes   : The structure containing all the nodes
% * parent      : The parent of the node to be checked [string]
% 
%% Outputs: 
% * checktot  : The confirmation whether the node is a valid child. This
%               should be 1 if it is and 0 otherwise 
% 
%% Author(s): Aram Vroom (2016)
% Email:  aram.vroom@strath.ac.uk

% Retreive the boundaries set for the node validity checks
CheckBounds = Inputs.NodeCheckBoundaries;

% Max dV_dep check
check1 = (norm(Attributes.dV_dep) <= CheckBounds(1));

% Minimum perihelion check
check2 = ((Attributes.kep_trans.a*(1-Attributes.kep_trans.e)) >= CheckBounds(2));

% Low-Thrust Check
check3 = ((Attributes.t_arr-ListNodes.(parent).attributes.t_arr -20)*86400*1e-7 >= CheckBounds(3)*(norm(Attributes.dV_dep)));

% Waiting time check
check4 = Attributes.t_dep-ListNodes.(parent).attributes.t_arr < CheckBounds(4);

% dV so far check
check5 = Attributes.dV_tot < CheckBounds(5);

% Confirm whether all conditions are met
checktot = check1*check2*check3*check4*check5;


end

