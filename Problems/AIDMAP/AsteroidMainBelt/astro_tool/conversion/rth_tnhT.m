function x_tnh = rth_tnhT(x_rth, a, e, f, mu)

% Vector reference frame Transformation.
% Radial-transversal-h reference frame to tangent-normal-h reference frame.
%
%	x_tnh = rth_tnhT(x_rth, a, e, f, mu)
%
% rth reference frame: {r, t, h}
%   r-axis: direction of the orbit radius
%   h-axis: direction of angular momentum
%   t-axis: in the orbit plane, completes the reference frame
% tnh reference frame: {t, n, h}
%   t-axis: tangent to the motion
%   h-axis: direction of angular momentum
%   n-axis: inward normal to t, in the orbit plane
%
% INPUT:
%        x_rth = vector to be transformed, expressed in {r, t, h}
%        a = semi-major axis
%        e = eccentricity
%        f = true anomaly from the pericentre [rad]
%        mu = gravitational constant
%
% OUTPUT:
%        x_tnh = vector transformed into {t, n, h}
%
% EXAMPLE:
%   Given a spacecraft in orbit:
%       - we have the thrust vector in {r, t, h};
%       - we want the thrust vector in {t, n, h}.
%   In this case:
%       x_rth = Thrust vector in {r, t, h};
%       a, e, f = Orbital parameters of the spacecraft;
%       mu = Gravitational constant of the attractor;
%       x_tnh = Thrust vector, transformed in {t, n, h}.
%
% FUNCTIONS CALLED: none
%
% - Camilla Colombo - 10/03/2006
%                   - 07/02/2006 - removed v from the inputs
% - Revised by Matteo Ceriotti - 10/01/2007
% - Matteo Ceriotti - 11/02/2008: Help improved.
%
% ------------------------- - SpaceART Toolbox - --------------------------

x_rth = x_rth(:);
p = a*(1-e^2);
n = sqrt(mu/a^3);
h = n*a^2*sqrt(1-e^2);
r = p/(1+e*cos(f));
v = sqrt(2*mu/r - mu/a);

sinb = h*e/(p*v)*sin(f);
cosb = h/(p*v)*(1+e*cos(f));

x_tnh = [sinb -cosb 0; cosb sinb 0; 0 0 1]'*x_rth;

return;