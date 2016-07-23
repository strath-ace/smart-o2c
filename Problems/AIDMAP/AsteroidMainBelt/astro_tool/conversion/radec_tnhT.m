function x_tnh = radec_tnhT(x_radec)

% Vector reference frame Transformation.
% Right Ascension and DEClination reference frame to tangentential-normal-h
% reference frame.
%
%	x_tnh = radec_tnhT(x_radec)
%
% radec reference frame: {r, alpha, delta}
%   right ascension and declination (spherical equatorial) reference frame:
%   r = modulus of the vector
%   alpha = in-plane right ascention angle, counted from the tangential
%           direction to the projection of the vector on the orbital plane
%           [rad]
%   delta = out-of-plane declination angle from the projection of the
%           vector on the orbital plane up to the vector itself [rad]
% tnh reference frame: {t, n, h}
%   t-axis: tangent to the motion
%   h-axis: direction of angular momentum
%   n-axis: inward normal to t, in the orbit plane
%
% INPUT:
%        x_radec = vector to be transformed, expressed in {r, alpha, delta}
%
% OUTPUT:
%        x_tnh = vector transformed into {t, n, h} (column)
%
% FUNCTIONS CALLED: none
%
% - Camilla Colombo - 07/12/2007
% - Revised by Matteo Ceriotti - 12/12/2007
% - Matteo Ceriotti - 11/02/2008: Help improved.
%
% ------------------------- - SpaceART Toolbox - --------------------------

x_tnh = [0;0;0];

x_tnh(1) = x_radec(1) * cos(x_radec(3))*cos(x_radec(2));
x_tnh(2) = x_radec(1) * cos(x_radec(3))*sin(x_radec(2));
x_tnh(3) = x_radec(1) * sin(x_radec(3));

return