function [Attributes] = MyAttributeCalcsMainBelt(Inputs, Parent, cityname, Attributes)
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%-----------Copyright (C) 2016 University of Strathclyde-------------
%
%
%
%% MyAttributeCalcsMainBelt: This function calculates the attributes that can not be retrieved from the unique ID.
% 
%% Inputs:
% * Inputs        : Structure containing the PhysarumSolver inputs
% * Parent        : Structure containing the parent
% * cityname      : The asteroid/planet of the node to be created [string]
% * Attributes    : The structure with attributes of the node to be created 
%                   known so far
% 
%% Outputs: 
% * Attributes   : The updated attributes [structure]
% 
%% Author(s): Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk

% Retrieve the data on the asteroids
Asteroids = Inputs.AdditionalInputs{1};

% Retrieve the cartesian coordiantes of the asteroid/planet at the arrival epoch
[Attributes.r_arr, Attributes.v_body] = StardustTool.CartesianElementsAt(Asteroids.(cityname), Attributes.t_arr);

% Save the body's orbital parameters
Attributes.kep_body = StardustTool.KeplerianElementsAt2(Asteroids.(cityname), Attributes.t_arr);

% Find the departure time by subtracting the ToF from the arrival time
Attributes.t_dep = Attributes.t_arr - Attributes.tof;

% Check if the node has a parent (excludes the root)
if ~isempty(Parent)

    % Add the ToF for this transfer to the total ToF of the parent
    % to find the new total ToF
    Attributes.tof_tot = Parent.attributes.tof_tot + Attributes.tof;
else

    % If no parent has been set, the total ToF is the current ToF
    Attributes.tof_tot = Attributes.tof;

    % The current keplerian orbit is then assumed to be that of the root
    kep = Asteroids.(cityname).getKeplerianElements;
    Attributes.kep_trans = CelestialBody('Transfer Orbit', kep.a, kep.e, kep.i, kep.OM, kep.W, kep.M0, kep.t0);
end 

end
