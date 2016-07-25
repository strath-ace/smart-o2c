% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%-----------Copyright (C) 2016 University of Strathclyde-------------
%


function [r, v] = EphSS(n, t)

% Ephemerides of the solar system.
% 
% 	[r, v] = EphSS(n, t)
% 
% It uses uplanet for planets, moon_eph for the Moon, and NeoEphemeris for
% asteroid ephemerides. Outputs cartesian position and velocity of the body, 
% centered in the Sun for all the bodies but the Moon (for which a
% cartesian Earth-centered reference frame is chosen).
% 
% INPUT
%  n   ID of the body:
%          1 to 9: Planets
%          11: Moon
%          >=12: NEOs
%  t   Time [d, MJD2000]. That is:
%      modified Julian day since 01/01/2000, 12:00
%      (MJD2000 = MJD-51544.5)
% 
% OUTPUT
%  r   Cartesian position of the body (Sun-centered for all bodies, 
%      Earth-centered for the Moon).
%  v   Cartesian velocity.
% 
% FUNCTIONS CALLED
%  uplanet, moon_eph, NeoEphemeris, kep2cart, astro_constants
% 
% - Matteo Ceriotti, 10-01-2007
%                    12-02-2007
% - Nicolas Croisard 03-05-2008 (moon_eph instead of uplanet for the Moon)
% - Revised by Camilla Colombo - 03-05-2008
% - Juan Manuel Romero Martin 10-12-2013 Modify the Astro_const to
%                                       AstroConstants Class. 
% 
% ------------------------- - SpaceART Toolbox - --------------------------

if n<11 % uplanet needed
    kep = uplanet(t, n);
elseif n==11 % moon_eph needed
    [r, v] = moon_eph(t); % Returns the cartesian position and velocity
else % NeoEphemeris needed
    kep = NeoEphemeris(t, n);
end

if n~=11 % Planet or asteroid, Sun-centered
    % Gravitational constant of the Sun
    % mu = astro_constants(4); % 132724487690 copied here for speed
    mu = AstroConstants.Sun_Planetary_Const; % [JMRM]
    
    % transform from Keplerian to cartesian
    car = kep2cart(kep, mu);
    r   = car(1:3);
    v   = car(4:6);
end

return