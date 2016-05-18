function x_tnh = car_tnhT(x_car,s_car)

% Vector reference frame Transformation.
% Cartesian reference frame to tangent-normal-h reference frame.
%
%   x_tnh = car_tnhT(x_car,s_car)
%
% car reference frame: {x,y,z}
%   inertial reference frame
% tnh reference frame: {t,n,h}
%   t-axis: tangent to the motion
%   h-axis: direction of angular momentum
%   n-axis: inward normal to t, in the orbit plane
%
% INPUT:
%        x_car = vector expressed in {x,y,z}
%        s_car = state vector [position, velocity] of the orbiting body,
%                in {x,y,z}
%
% OUTPUT:
%        x_tnh = vector transformed into {t,n,h}
%
% EXAMPLE:
%   Given a spacecraft in orbit:
%       - we have the thrust vector in {x,y,z};
%       - we want the thrust vector in {r,n,h}.
%   In this case:
%       x_car = Thrust vector in {x,y,z};
%       s_car = [position, velocity] of the spacecraft in {x,y,z};
%       x_tnh = Thrust vector, transformed into {t,n,h}.
%
% FUNCTIONS CALLED: none
%
% - Camilla Colombo - 03/03/2006
% - Revised by Matteo Ceriotti - 10/01/2007
% - Matteo Ceriotti - 11/02/2008: Help improved.
%
% ------------------------- - SpaceART Toolbox - --------------------------

x_car = x_car(:);
s = s_car(:);
r = s(1:3);
v = s(4:6);
t_ = v/norm(v);
h = cross(r,v);
h_ = h/norm(h);
n_ = cross(h_,t_);

A = [t_ n_ h_];

x_tnh = A'*x_car;

return