function v1=euler_axis_angle(v,n,theta)

% Rotates a vector about an axis of a given angle (Eules axis and angle
% rotation).
%
% v1 = euler_axis_angle(v, n, theta)
%
% Note: it uses the right-hand rule. (See Analytical Mechanics of Space
%       Systems, Schaub< Junkins, pp. 90)
%
% INPUT
%   v       Vector to be rotated.
%   n       Axis of rotation.
%   theta   Angle of rotation [rad].
%
% OUTPUT
%   v1      Rotated vector.
%
% FUNCTIONS CALLED
%   (none)
%
% Matteo Ceriotti, 11-01-2007
% Revised by Camilla Colombo - 14-05-2007
%
% ------------------------- - SpaceART Toolbox - --------------------------

v=v(:);
n=n/norm(n);
R=[cos(theta)+(1-cos(theta))*n(1)^2, (1-cos(theta))*n(1)*n(2)+sin(theta)*n(3), (1-cos(theta))*n(1)*n(3)-sin(theta)*n(2);...
    (1-cos(theta))*n(1)*n(2)-sin(theta)*n(3), cos(theta)+(1-cos(theta))*n(2)^2, (1-cos(theta))*n(2)*n(3)+sin(theta)*n(1);...
    (1-cos(theta))*n(1)*n(3)+sin(theta)*n(2), (1-cos(theta))*n(2)*n(3)-sin(theta)*n(1), cos(theta)+(1-cos(theta))*n(3)^2];
v1=R*v;