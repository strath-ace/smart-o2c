function [] = InitialiseAsteroidsMainBelt(epoch_start, epoch_end, filenames)
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%-----------Copyright (C) 2018 University of Strathclyde and Authors-----------
%
%
%
%% InitialiseAsteroidsMainBelt: This function retrieves the asteroid names and the passing epochs through the nodal points required for the Main Belt problem
% 
%% Inputs:
% * epoch_start : The epoch from which the nodal points are to be
%                 calculated [MJD2000 date]
% * epoch_end   : The epoch up to which the nodal points should be
%                 calculated [MJD2000 date]
% * filenames   : Structure that contains the file names
% 
%% Outputs: 
% * 
% 
%% Author(s): Aram Vroom (2016)
% Email:  aram.vroom@strath.ac.uk

if isunix
    addpath(genpath(strcat(fileparts(pwd),'/astro_tool')));
else
    addpath(genpath(strcat(fileparts(pwd),'\astro_tool')));
end

% Obtain the Asteroids structure that contains the CelestialBody objects for
% each asteroid
[Asteroids] = EvaluateAsteroidsXLS(filenames.AsteroidsFileName, filenames.MatFileName, filenames.NameFile);

% Retrieve their names
asteroidnames = fieldnames(Asteroids);

% Loop over all the asteroids
for i = 1:length(asteroidnames)
    
    % Ensure the name is written as characters
    asteroid = char(asteroidnames(i));

    % Retrieve the gravitational parameter of the sun
    mu = AstroConstants.Sun_Planetary_Const;

    % Generate a CelestialBody that holds the virtual orbit      
    virtorbit = CelestialBody('virtualorbit', ...        % Name
                                1, ...                   % a   Semimajor axis [AU]   
                                0.035999947927132, ...   % e   Eccentricity 
                                0, ...                   % i   Inclination [deg]
                                -11.26064, ...           % OM  Asc. Node/raan [deg]
                                102.94719, ...           % W   Arg. Perigee [deg]
                                0, ...                   % M0  Mean anomoly, M at time given t0 [deg]
                                0);                     % t0  Time at which Mo is given [MJD2000]  

    % Retrieve the asteroid's orbit
    asteroidorbit = Asteroids.(asteroid).getKeplerianElements();

    % Compute Nodal points
    [Mnode1(i), Mnode2(i), ~, ~, ~, ~, ~] = computeNodalPoints_M0(virtorbit, Asteroids.(asteroid), mu);
    
    % Compute the passing epochs
    epochsnode1 = computeNodalPassesEpochs(asteroidorbit, Mnode1(i), epoch_start, epoch_end, mu);
    epochsnode2 = computeNodalPassesEpochs(asteroidorbit, Mnode2(i), epoch_start, epoch_end, mu);

    % Save the passing epochs into a single array
    epochsnode{i, 1} = [epochsnode1 epochsnode2];

end

% Save the epochs
save(filenames.epochsnodename, 'epochsnode')
end



