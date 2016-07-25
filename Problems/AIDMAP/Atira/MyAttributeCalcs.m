% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%-----------Copyright (C) 2016 University of Strathclyde-------------
%


function [Attributes] = MyAttributeCalcs(Inputs, Parent, cityname, Attributes)
%% MyAttributeCalcs: This function calculates the attributes that can not be retrieved from the unique ID.
% 
%% Inputs:
% * Inputs        : Structure containing the inputs defined by the
%                   user. This variable is currently not used in this file, 
%                   but shown here to illustrate the ability for users 
%                   to use this structure.
% * Parent        : Structure containing the parent
% * cityname      : The asteroid/planet of the node to be created [string]
% * Attributes    : The structure with attributes of the node to be created 
%                   known so far
% 
%% Outputs: 
% * Attributes   : The updated attributes structure
% 
%% Author(s): Marilena Di Carlo (2014), Aram Vroom (2016)
% Email:  marilena.di-carlo@strath.ac.uk, aram.vroom@strath.ac.uk

% Check if the asteroid/planet is in the Asteroid class
if ismember(cityname, fieldnames(Asteroids))
    
    % Retrieve the cartesian coordiantes of the asteroid/planet at the arrival epoch
    Attributes.r_arr = StardustTool.CartesianElementsAt(Asteroids.(cityname), Attributes.t_arr);
        
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

% Check if the asteroid/planet is a planet
elseif ismember(cityname, fieldnames(Planets));
    
    % Obtain the arrival coordinates from the planetary ephemeris
    [Attributes.r_arr, ~] = EphSS(Planets.(cityname), Attributes.t_arr);
    
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
        
        kep_elements = EphSS_kep(Planets.(cityname), Attributes.t_dep);              
        
        % True anomaly at departure position over transfer orbit [rad]
        theta = kep_elements(6);
        ecc = kep_elements(2);

        % Compute the eccentric anomaly at the departure position [rad]
        cos_E = ( cos(theta) + ecc ) / (1 + ecc * cos(theta) );
        sin_E = ( sin(theta) * sqrt(1 - ecc^2) ) / (1 + ecc * cos(theta) );

        E = atan2(sin_E, cos_E);
        E = mod(E, 2*pi);

        % Compute the mean anomaly at the departure position [rad]
        M = E - ecc*sin(E); 

        % Mean anomaly at the departure position [deg]
        M = M * 180/pi;

        % Create kep_trans object
        au2km = AstroConstants.Astronomical_Unit;
        Attributes.kep_trans = CelestialBody('Transfer Orbit',       ... % Name of the CelestialBody (This case current trajectory
                                             kep_elements(1)/au2km,  ... % Semimajor axis [AU]
                                             kep_elements(2),        ... % Eccentricity
                                             kep_elements(3)*180/pi, ... % Inclination [deg]
                                             kep_elements(4)*180/pi, ... % Asc. Node/Raan [deg]
                                             kep_elements(5)*180/pi, ... % Arg. Perigee [deg]
                                             M,                      ... % Mean anomoly, M at time given t0 [deg]
                                             Attributes.t_dep);          % Time at which Mo is given [MJD2000]  

   end 
    
else
    
    % If the asteroid/planet is not in the Asteroids and Planets list, set the r_arr, tof_tot and
    % t_dep to default values.
    Attributes.r_arr = zeros(1, 3);
    Attributes.tof_tot = 0;
    Attributes.t_dep = 0;
end


end
