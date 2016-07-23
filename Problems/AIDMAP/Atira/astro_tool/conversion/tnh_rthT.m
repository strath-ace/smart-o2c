function x_rth = tnh_rthT(x_tnh, a, e, f, mu)

% Vector reference frame Transformation.
% Tangent-normal-h reference frame to radial-transversal-h reference frame.
%
%   x_rth = tnh_rthT(x_tnh, a, e, f, mu)
%
% tnh reference frame: {t, n, h}
%   t-axis: tangent to the motion
%   h-axis: direction of angular momentum
%   n-axis: inward normal to t, in the orbit plane
% rth reference frame: {r, t, h}
%   r-axis: direction of the orbit radius
%   h-axis: direction of angular momentum
%   t-axis: in the orbit plane, completes the reference frame
%
% INPUT:
%        x_tnh = vector to be transformed, expressed in {t, n, h}
%        a = semi-major axis
%        e = eccentricity
%        f = true anomaly from the perigee [rad]
%        mu = gravitational constant
%
% OUTPUT:
%        x_rth = vector transformed into {r, t, h}
%
% EXAMPLE:
%   Given a spacecraft in orbit:
%       - we have the thrust vector in {t, n, h};
%       - we want the thrust vector in {r, t, h}.
%   In this case:
%       x_tnh = Thrust vector in {t, n, h};
%       a, e, f = Orbital parameters of the spacecraft;
%       mu = Gravitational constant of the attractor;
%       x_rth = Thrust vector, transformed in {r, t, h}.
%
% FUNCTIONS CALLED: none
%
% - Camilla Colombo - 23/02/2006
%                   - 07/02/2006 - removed v from the inputs
% - Revised by Matteo Ceriotti - 10/01/2007
% - Matteo Ceriotti - 11/02/2008: Help improved.
%
% ------------------------- - SpaceART Toolbox - --------------------------

x_tnh = x_tnh(:);
p = a*(1-e^2);
n = sqrt(mu/a^3);
h = n*a^2*sqrt(1-e^2);
r = p/(1+e*cos(f));
v = sqrt(2*mu/r - mu/a);

sinb = h*e/(p*v)*sin(f);
cosb = h/(p*v)*(1+e*cos(f));

x_rth = [sinb -cosb 0; cosb sinb 0; 0 0 1]*x_tnh;

return