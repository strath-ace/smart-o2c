% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%-----------Copyright (C) 2016 University of Strathclyde-------------
%


function [epochs, error_flag] = computeNodalPassesEpochs( keporb, M, epoch_start, epoch_end, mu)
%% computeNodalPassesEpochs: computes the passing epochs at the given Mean Anomaly M0 within the given time domain [epoch_start, epoch_end].
% 
%% Inputs:
%  keporb       : Keplerian Elements for the current Orbit. The keporb has
%                 to be a KeplerianElements Object that contains the 
%                 following attributes:
% 
%                   a   Semimajor axis [AU]     
%                   e   Eccentricity
%                   i   Inclination [deg]     
%                   OM  Asc. Node/raan [deg]
%                   W   Arg. Perigee [deg]     
%                   M0  Mean anomoly, M at time given t0 [deg]
%                   t0  Time at which Mo is given [MJD2000]                 
%                 
%  M            : Mean Anomaly [deg]
%  epoch_start  : Start Epoch for the Time Frame [MDJ2000]
%  epoch_end    : End Epoch for the Time Frame [MDJ2000]
%  mu           : Gravitational Parameter [km^3/s^2]
% 
%% Outputs:
% epochs        : The passing epochs
% error_flag    : Flag that indicates whether an error occured
% 
%% Author(s): Juan Manuel Romero Martin (2014)
% Email:  juan.romero-martin@strath.ac.uk
% 
% Equation used to compute the Epoch at given position:
% 
% M = n * (t - t0) + M0
% |
% |--> t = (M - M0)/n + t0;
% 
% where: 
%        n = sqrt(mu/a^3)
% 
% --------------------------------------------------------------------------

% Set the error flag to true
error_flag = 1;

% Sanity Check: The keporb is given in Keplerian Object
if isa(keporb, 'KeplerianElements') == 0    
    error('The keporb must be a KeplerianElements Object.')
end

% Sanity Check: First Obirt is given in Keplerian Object
if epoch_start > epoch_end 
    error('The Start Epoch can NOT be earlier than End Epcoh.')
end

% Sanity Check: Valid Mean Anomaly
if isnan(M)
    error('The Mean Anomaly is Not A Number  (NaN).')
end

% Convert the Semimajor from AU to KM
a = keporb.a * AstroConstants.Astronomical_Unit;

% Compute the mean motion [1/sec]
n = sqrt(mu / a^3);

% Compute the period [sec]
T = 2*pi/n;

% Convert the period in days [day]
T_days = T/60/60/24;

% Compute the first passing Epoch for the given Mean Anomaly position [MDJ2000]
epoch_first = (M*pi/180 - keporb.M0*pi/180)/(n*86400) + keporb.t0;

% Sanity Check: check if the Start Epoch is later than the computed first 
% epoch
if epoch_start >= epoch_first
    
    % Compute the time difference from the first passing epoch
    epoch_diff_start = epoch_start - epoch_first;
    
    % Get the number of revolution need to reach the start epoch
    epoch_first_start = epoch_first + ceil(epoch_diff_start / T_days) * T_days;
        
end

% Sanity Check: check if the Start Epoch is too earlier than the computed 
% first epoch
if epoch_first >= epoch_start
    
    % Compute the time difference from the first passing epoch
    epoch_diff_start = epoch_first - epoch_start;
    
    % Get the number of revolution need to reach the start epoch
    epoch_first_start = epoch_first - floor(epoch_diff_start / T_days) * T_days;
    
end

% Compute the number of possible passing epochs during the given time
% domain
num_rev = floor(epoch_end - epoch_first_start)/T_days;

% Pre-allocate memory to save the passing epochs available for the given 
% time domain 
epochs = 1:num_rev;

% Compute the passing epochs
for i = 0 : num_rev    
    epochs(i + 1) = epoch_first_start + T_days*i;    
end

% Set the error flag to false
error_flag = 0;