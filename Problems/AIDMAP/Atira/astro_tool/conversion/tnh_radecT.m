function x_radec = tnh_radecT(x_tnh)

% Vector reference frame Transformation.
% Tangent-normal-h reference frame to Right Ascension and DEClination
% reference frame.
%
%   x_radec = tnh_radecT(x_tnh)
%
% tnh reference frame: {t, n, h}
%   t-axis: tangent to the motion
%   h-axis: direction of angular momentum
%   n-axis: inward normal to t, in the orbit plane
% radec reference frame: {r, alpha, delta}
%   right ascension and declination (spherical equatorial) reference frame:
%   r = modulus of the vector
%   alpha = in-plane right ascention angle, counted from the tangential
%           direction to the projection of the vector on the orbital plane
%           [rad]
%   delta = out-of-plane declination angle from the projection of the
%           vector on the orbital plane up to the vector itself [rad]
%
% INPUT:
%        x_tnh = vector to be transformed, expressed in {t, n, h}
%
% OUTPUT:
%        x_radec = vector transformed into {r, alpha, delta} (column)
%
% FUNCTIONS CALLED: none
%
% - Camilla Colombo - 23/10/2007
% - Revised by Matteo Ceriotti - 12/12/2007
% - Matteo Ceriotti - 11/02/2008: Help improved.
%
% ------------------------- - SpaceART Toolbox - --------------------------

x_radec = [0;0;0];

x_radec(1) = (x_tnh(1)^2+x_tnh(2)^2+x_tnh(3)^2)^0.5;
x_radec(2) = atan2(x_tnh(2), x_tnh(1));
x_radec(3) = asin(x_tnh(3)/x_radec(1));

return