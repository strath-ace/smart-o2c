function kep = EphSS_kep(n, t)
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
<<<<<<< HEAD
%-----------Copyright (C) 2018 University of Strathclyde and Authors-----------
=======
%-----------Copyright (C) 2016 University of Strathclyde-------------
>>>>>>> 5b7361d93c9119cf1d2e9e6c885bed93f924d71b
%
%
%
%% EphSS_kep: This file obtains the ephemerides of the solar system in Keplerian parameters.
% 
%  kep = EphSS_kep(n, t)
% 
% Uses uplanet for planets, moon_eph for the Moon, and NeoEphemeris for
% asteroid ephemerides. Outputs the keplerian elements of the body, 
% centered in the Sun for all the bodies but the Moon (Earth-centered).
% 
%% Inputs:
%  n   ID of the body:
%          1 to 9: Planets
%          11: Moon
%          >=12: NEOs
%  t   Time [d, MJD2000]. That is:
%      modified Julian day since 01/01/2000, 12:00
%      (MJD2000 = MJD-51544.5)
% 
%% Outputs:
%  kep     Keplerian elements of the body (Sun-centered for all bodies, 
%          Earth-centered for the Moon).
%          kep = [a e i Om om theta] [km, rad]
% 
% FUNCTIONS CALLED
%  uplanet, moon_eph, NeoEphemeris, cart2kep, astro_constants
% 
% - Nicolas Croisard 03-05-2008
% - Revised by Matteo Ceriotti - 03-05-2008
% 
% ------------------------- - SpaceART Toolbox - --------------------------

if n<11 % uplanet needed
    kep = uplanet(t, n);
elseif n==11 % moon_eph needed
    [r, v] = moon_eph(t); % Returns the cartesian position and velocity
    mu = astro_constants(13); % Gravitational constant of the Earth
    kep = cart2kep([r, v], mu); % Transform from cartesian to Keplerian 
else % NeoEphemeris needed
    kep = NeoEphemeris(t, n);
end

return
