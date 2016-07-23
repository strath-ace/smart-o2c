function x_car = rth_carT(x_rth, s_car)

% Vector reference frame Transformation.
% Radial-trasversal-h reference frame to Cartesian reference frame.
%
%	x_car = rth_carT(x_rth, s_car)
%
% rth reference frame: {r, t, h}
%   r-axis: direction of the orbit radius
%   h-axis: direction of angular momentum
%   t-axis: in the orbit plane, completes the reference frame
% car reference frame: {x, y, z}
%   inertial reference frame
%
% INPUT:
%        x_rth = vector to be transformed, expressed in {r, t, h}
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

x_rth = x_rth(:);
s = s_car(:);
r = s(1:3);
v = s(4:6);
r_ = r/norm(r);
h = cross(r, v);
h_ = h/norm(h);
t_ = cross(h_, r_);

A = [r_ t_ h_];

x_car = A*x_rth;

return