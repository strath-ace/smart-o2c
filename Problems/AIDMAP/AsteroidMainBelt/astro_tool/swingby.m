function [v2,delta]=swingby(v1,rp,mu,gamma,n_r)
% Computes the outgoing velocity after the swingby of a planet, given
% incoming velocity and swingby parameters.
%
% v2 = swingby(v1, rp, mu, gamma, n_r);
%
% See the documentation for a detailed explaination of the angles.
%
% INPUT
%   v1      Incoming velocity, before the swingby, relative to the planet.
%   rp      Radius of pericentre of the hyperbola.
%   mu      Gravitational constant of the planet.
%   gamma   Plane angle [rad]. This angle identifies the inclination of the
%           hyperbola plane around the incoming velocity vector, and is the
%           angle between the vector n_r and the vector normal to the
%           hyperbola plane.
%   n_r     Reference vector, used as a origin to measure gamma. In
%           principle, this vector is arbitrary.
%           A choice can be to use the normal to the plane containing the
%           incoming velocity and the heliocentric velocity of the planet:
%               n_r = cross(v1, vp) / norm(cross(v1, vp));
%           being vp the heliocentric velocity of the planet.
%
% OUTPUT
%   v2      Outgoing velocity, after the swingby, relative to the planet.
%
% FUNCTION CALLED
%   euler_axis_angle
%
% Matteo Ceriotti, 11-01-2007
% Modified by: Matteo Ceriotti, 14-05-2007: Help improved.
% Revised by: Camilla Colombo, 15-05-2007
%
% ------------------------- - SpaceART Toolbox - --------------------------

mu_rp=mu/rp;
theta_inf=acos(-mu_rp/(norm(v1)^2+mu_rp)); % See Kaplan "Modern spacecraft dynamics and control", pag. 93
delta=2*theta_inf-pi; % Deflection angle

n_pi=euler_axis_angle(n_r,v1,gamma); % Rotates n_r around v1
v2=euler_axis_angle(v1,n_pi,delta); % Rotates v1 around n_pi
