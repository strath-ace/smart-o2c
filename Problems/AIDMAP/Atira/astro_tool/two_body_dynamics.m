function y=two_body_dynamics(t, s, mu)
% two_body_dynamics - Dynamics of the two-body problem
%
%   y = two_body_dynamics(t, s, mu)
%
% Computes the derivative of the state vector (position and velocity for
% the two-body problem dynamics.
%
% INPUT
%   t   Time
%   s   State vector (position and velocity)
%   mu  Planetary constant.
%
% OUTPUT
%   y   Derivative of the state vector
%
% EXAMPLE
%   One Earth orbit:
%       [x0, v0]=EphSSfun(3, 0);
%       [t, s]=ode45(@two_body_dynamics, [0, 365*86400], [x0, v0], [], astro_constants(4));
%       plot3(s(:, 1), s(:, 2), s(:, 3));
%
% Matteo Ceriotti, 2006
% Revised by Camilla Colombo - 06-02-2007
%
% ------------------------- - SpaceART Toolbox - --------------------------

y = [       s(4:6)       ;...
     -mu/norm(s(1:3))^3*s(1:3)];

return