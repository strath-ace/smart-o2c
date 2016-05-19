function x_radec = car_radecT(x_car,s_car)

% Vector reference frame Transformation.
% Cartesian reference frame to Right Ascension and DEClination reference
% frame.
%
%   x_radec = car_radecT(x_car,s_car)
%
% car reference frame: {x,y,z}
%   inertial reference frame
% radec reference frame: {r,alpha,delta}
%   right ascension and declination (spherical equatorial) reference frame:
%   r = modulus of the vector
%   alpha = in-plane right ascention angle, counted from the tangential
%           direction to the projection of the vector on the orbital plane
%           [rad]
%   delta = out-of-plane declination angle from the projection of the
%           vector on the orbital plane up to the vector itself [rad]
%
% INPUT:
%        x_car = vector to be transformed, expressed in {x,y,z}
%        s_car = state vector [position, velocity] of the orbiting body,
%                in {x,y,z}
%
% OUTPUT:
%        x_radec = vector transformed into {r,alpha,delta}
%
% EXAMPLE:
%   Given a spacecraft in orbit:
%       - we have the thrust vector in {x,y,z};
%       - we want the thrust vector in {r,alpha,delta}.
%   In this case:
%       x_car = Thrust vector in {x,y,z};
%       s_car = [position, velocity] of the spacecraft in {x,y,z};
%       x_radec = Thrust vector, transformed into {r,alpha,delta}.
%
% FUNCTIONS CALLED: car_tnhT,tnh_radecT
%
% - Camilla Colombo - 23/10/2007
% - Revised by Matteo Ceriotti - 12/12/2007
% - Matteo Ceriotti - 11/02/2008: Help improved.
%
% ------------------------- - SpaceART Toolbox - --------------------------

x_tnh = car_tnhT(x_car,s_car);
x_radec = tnh_radecT(x_tnh);

return