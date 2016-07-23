function x_car = radec_carT(x_radec, s_car)

% Vector reference frame Transformation.
% Right Ascension and DEClination reference frame to Cartesian reference
% frame.
%
%   x_car = radec_carT(x_radec, s_car)
%
% radec reference frame: {r, alpha, delta}
%   right ascension and declination (spherical equatorial) reference frame:
%   r = modulus of the vector
%   alpha = in-plane right ascention angle, counted from the tangential
%           direction to the projection of the vector on the orbital plane
%           [rad]
%   delta = out-of-plane declination angle from the projection of the
%           vector on the orbital plane up to the vector itself [rad]
% car reference frame: {x, y, z}
%   inertial reference frame
%
% INPUT:
%        x_radec = vector to be transformed, expressed in {r, alpha, delta}
%        s_car = state vector [position, velocity] expressed in {x, y, z}
%
% OUTPUT:
%        x_car = vector transformed into {x, y, z}
%
% EXAMPLE:
%   Given a spacecraft in orbit:
%       - we have the thrust vector in {r, alpha, delta};
%       - we want the thrust vector in {x, y, z}.
%   In this case:
%       x_radec = Thrust vector in {r, alpha, delta};
%       s_car = [position, velocity] of the spacecraft in {x, y, z};
%       x_car = Thrust vector, transformed in {x, y, z}.
%
% FUNCTIONS CALLED: radec_tnhT, tnh_carT
%
% - Camilla Colombo - 07/12/2007
% - Revised by Matteo Ceriotti - 12/12/2007
% - Matteo Ceriotti - 11/02/2008: Help improved.
%
% ------------------------- - SpaceART Toolbox - --------------------------

x_tnh = radec_tnhT(x_radec);
x_car = tnh_carT(x_tnh, s_car);

return