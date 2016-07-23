function x_rth = car_rthT(x_car, s_car)

% Vector reference frame Transformation.
% Cartesian reference frame to radial-trasversal-h reference frame.
%
%	x_rth = car_rthT(x_car, s_car)
%
% car reference frame: {x, y, z}
%   inertial reference frame
% rth reference frame: {r, t, h}
%   r-axis: direction of the orbit radius
%   h-axis: direction of angular momentum
%   t-axis: in the orbit plane, completes the reference frame (inward)
%
% INPUT:
%        x_car = vector to be transformed, expressed in {x, y, z}
%        s_car = state vector [position, velocity] of the orbiting body, 
%                in {x, y, z}
%
% OUTPUT:
%        x_rth = vector transformed into {r, t, h}
%
% EXAMPLE:
%   Given a spacecraft in orbit:
%       - we have the thrust vector in {x, y, z};
%       - we want the thrust vector in {r, t, h}.
%   In this case:
%       x_car = Thrust vector in {x, y, z};
%       s_car = [position, velocity] of the spacecraft in {x, y, z};
%       x_rth = Thrust vector, transformed into {r, t, h}.
%
% FUNCTIONS CALLED: none
%
% - Camilla Colombo - 03/03/2006
% - Revised by Matteo Ceriotti - 10/01/2007
% - Matteo Ceriotti, Nicolas Croisard - 24/01/2008: Help improved.
%
% ------------------------- - SpaceART Toolbox - --------------------------

x_car = x_car(:);
s = s_car(:);
r = s(1:3);
v = s(4:6);
r_ = r/norm(r);
h = cross(r, v);
h_ = h/norm(h);
t_ = cross(h_, r_);

A = [r_ t_ h_];

x_rth = A'*x_car;

return