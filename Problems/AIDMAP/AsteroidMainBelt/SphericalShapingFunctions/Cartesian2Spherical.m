function [pos,vel] = Cartesian2Spherical(r,v)
% Conversion from cartesian coordiantes (x,y,z and velocities) to spherical
% coordinates (r,theta,phi and velocities)
% Note: r,v vectors must be column vectors; not working if r,v are matrix
% Note: use atan2 for theta to find the right quadrant!
% Angles in radians

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Niccolo' Gastaldello, October 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute pos,vel

rad = sqrt(r(1)^2+r(2)^2+r(3)^2);
theta = atan2(r(2),r(1));
phi = asin(r(3)/rad);

pos = [rad;theta;phi];

vel = [cos(phi)*cos(theta)    cos(phi)*sin(theta)    sin(phi);
           -sin(theta)             cos(theta)             0;
      -sin(phi)*cos(theta)    -sin(phi)*sin(theta)  cos(phi)]*v;

end

