function yend = plot_trajectory(r0,v0,tof,mu,s)
% Plots a keplerian orbit of a body around a cental mass, given the initial
% position and velocity in cartesian coordinates, and the final time.
%
%   plot_trajectory(r0, v0, tof, mu, s)
%
% Units of measure consistent each other.
%
% INPUT
%   r0      Cartesian initial position.
%   v0      Cartesian initial velocity.
%   tof     Time of flight.
%   mu      Gravity constant of the central body.
%   s       Character string made from one element from the following
%           column to define the color of the line plotted
%               b     blue
%               g     green
%               r     red
%               c     cyan
%               m     magenta
%               y     yellow
%               k     black
%
% OUTPUT
%   (none)
%
% FUNCTIONS CALLED
%   two_body_dynamics
%
% Matteo Ceriotti, 06-02-2007
%
% Revised by Camilla Colombo - 06-02-2007
% Camilla Colombo, 24-11-2007 - s added as an input, yend added as output
%
% ------------------------- - SpaceART Toolbox - --------------------------

if nargin < 5
    s = 'b';
end

% Initialises figure
line(0,0,0,'Marker','*','Color','k');

% Initialise integrator used to plot trajectories
ode45options = odeset('abstol',1e-11,'reltol',1e-9);
[t,y]=ode45(@two_body_dynamics,[0, tof],[r0,v0],ode45options,mu);
yend = y(end,:);

line(y(:,1),y(:,2),y(:,3),'Color',s);
% line(y(1,1),y(1,2),y(1,3),'Marker','o','Color',s) % Marker at departure point
% line(yend(1),yend(2),yend(3),'Marker','o','Color',s) % Marker at arrival point
%plot3(y(:,1),y(:,2),y(:,3),s)
