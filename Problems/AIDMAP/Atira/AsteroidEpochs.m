% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%-----------Copyright (C) 2016 University of Strathclyde-------------
%


%% AsteroidEpochs: This script calculates the passing epochs of all the asteroids through their nodal points, given a time-frame. 
% It should be run first when a new set of Asteroids is used to obtain their nodal epochs.
% 
%% Inputs:
% * epoch_start : The epoch from which the nodal points are to be
%                 calculated [MJD2000 date]
% * epoch_end   : The epoch up to which the nodal points should be
%                 calculated [MJD2000 date]
% * Asteroids   : The file that contains the class which describes the
%                 asteroids
% 
%% Outputs: 
% * epochsnode  : Epochs at which the asteroid passes its nodal points
% 
%% Author(s): Aram Vroom (2016)
% Email:  aram.vroom@strath.ac.uk

% Initialise
close all; clear all; clc

%Add the astro_tool directory
if isunix
     addpath(genpath(strcat(fileparts(pwd), '/common')));
else
     addpath(genpath(strcat(fileparts(pwd), '\common')));
end

% Define input & output directory
SaveDir = 'IO_Dir\';

% If the user is using a Linux or Mac version of MATLAB, replace the
% backslashes by forward slashes in the Input / Output directory
if isunix
    SaveDir = strrep(SaveDir,'\','/');     
end

% Add this path
addpath(SaveDir);

% Input start and end epochs [MJD2000]
epoch_start = 7304.5;
epoch_end = 10957.5;

% Load the file containing the asteroids
Asteroids = Asteroids();

% Retrieve the Sun's standard gravitational parameter
mu = AstroConstants.Sun_Planetary_Const;

% Obtain the asteroid names
asteroidnames = fieldnames(Asteroids);

% Generate a CelestialBody that holds the virtual orbit        
virtorbit = CelestialBody('virtualorbit', ...        % Name
                            1, ...                   % a   Semimajor axis [AU]   
                            0.035999947927132, ...   % e   Eccentricity 
                            0, ...                   % i   Inclination [deg]
                            -11.26064, ...           % OM  Asc. Node/raan [deg]
                            102.94719, ...           % W   Arg. Perigee [deg]
                            0, ...                   % M0  Mean anomoly, M at time given t0 [deg]
                            0);                     % t0  Time at which Mo is given [MJD2000]  

% Loop over all the asteroids
for i = 1:length(asteroidnames)
    
% Convert the name of the asteorid currently being evaluated to characters
asteroid = char(asteroidnames(i));

% Retrieve the Keplerian elements
asteroidorbit = Asteroids.(asteroid).getKeplerianElements();

% Compute the location of Nodal points (the intersection between the two
% orbits)
[Mnode1, Mnode2, error_status, theta1, E1, theta2, E2 ] = computeNodalPoints_M0(virtorbit, Asteroids.(asteroid), mu);

% Find the passing epochs
epochsnode1 = computeNodalPassesEpochs(asteroidorbit, Mnode1, epoch_start, epoch_end, mu);
epochsnode2 = computeNodalPassesEpochs(asteroidorbit, Mnode2, epoch_start, epoch_end, mu);

% Save the passing epochs into a single array
epochsnode{i, 1} = sort([epochsnode1 epochsnode2]);

end

% Output the passing epoch into a file
save(strcat(SaveDir, 'epochsnode'), 'epochsnode')


