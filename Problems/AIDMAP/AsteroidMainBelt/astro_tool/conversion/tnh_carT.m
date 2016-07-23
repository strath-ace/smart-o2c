function x_car = tnh_carT(x_tnh, s_car)

% Vector reference frame Transformation.
% Tangent-normal-h reference frame to cartesian reference frame.
%
%	x_car = tnh_carT(x_tnh, s_car)
%
% tnh reference frame: {t, n, h}
%   t-axis: tangent to the motion
%   h-axis: direction of angular momentum
%   n-axis: inward normal to t, in the orbit plane
% car reference frame: {x, y, z}
%   inertial reference frame
%
% INPUT:
%        x_tnh = vector to be transformed, expressed in {t, n, h}
%        s_car = state vector [position, velocity] of the orbiting body, 
%                in {x, y, z}
%
% OUTPUT:
%        x_car = vector transformed into {x, y, z}
%
% EXAMPLE:
%   Given a spacecraft in orbit:
%       - we have the thrust vector in {r, t, h};
%       - we want the thrust vector in {x, y, z}.
%   In this case:
%       x_rth = Thrust vector in {r, t, h};
%       s_car = [position, velocity] of the spacecraft in {x, y, z};
%       x_car = Thrust vector, transformed in {x, y, z}.
%
% FUNCTIONS CALLED: none
%
% - Camilla Colombo - 03/03/2006
% - Revised by Matteo Ceriotti - 10/01/2007
% - Matteo Ceriotti - 11/02/2008: Help improved.
%
% ------------------------- - SpaceART Toolbox - --------------------------

x_tnh = x_tnh(:);
s = s_car(:);
r = s(1:3);
v = s(4:6);
t_ = v/norm(v);
h = cross(r, v);
h_ = h/norm(h);
n_ = cross(h_, t_);

A = [t_ n_ h_];

x_car = A*x_tnh;

return