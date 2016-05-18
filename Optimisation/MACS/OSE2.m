close all
clear all
clc

% Environmental parameters

rho = 1e-14;
grav_const = 398600*10^9;                                                   %[m3/s2]
earth_rad = 6371e3;                                                         %[m]

% Simulation parameters

tmax = 3600*2;                                                                % [s]

% All positions and velocity are ABSOLUTE (wrt fixed earth)

% Initial condition of chaser and its physical parameters
chaser.altitude = 90e3;                                                     %[m, from surface]
chaser.latitude = 0*pi/180;                                                 %[rad, from equator]
chaser.longitude = 0*pi/180;                                                %[rad, from Greenwich]

r = chaser.altitude+earth_rad;
x = r*cos(chaser.latitude)*cos(chaser.longitude);
y = r*cos(chaser.latitude)*sin(chaser.longitude);
z = r*sin(chaser.latitude);

chaser.x0 = [x;y;z];                                                        %[m]
chaser.xdot0 = [0;(grav_const/r^3)^0.5*r;0];                                %[m/s]
chaser.mass = 1435;                                                         %[kg]
chaser.inertia = [2040 130 25; 130 1670 -55; 25 -55 2570];                  %[kg/m2]
chaser.bbox = [3;5;2];                                                      %[m]
chaser.cd = [2.9;2.9;2.9];                                                  %[]
chaser.cm = [0;0;0];                                                        %[];
chaser.area = [56.64;56.64;56.64];                                           %[m^2];

% Initial condition of targer and its physical parameters
target.altitude = 200e3;                                                    %[m, from surface]
target.latitude = 0*pi/180;                                                 %[rad, from equator]
target.longitude = 0*pi/180;                                                %[rad, from Greenwich]

r = target.altitude+earth_rad;
x = r*cos(target.latitude)*cos(target.longitude);
y = r*cos(target.latitude)*sin(target.longitude);
z = r*sin(target.latitude);

target.x0 = [x;y;z];                                                        %[m]
target.xdot0 = [0;(grav_const/r^3)^0.5*r;0];                                %[m/s]
target.mass = 7792;                                                         %[kg]
target.inertia = [17012 401 -2167; 401 124725 345; -2167 345 129009];       %[kg/m2]
target.bbox = [4;5;27];                                                     %[m]
target.cd = [2.9;2.9;2.9];                                                  %[]
target.cm = [0;0;0];                                                        %[];
target.area = [56.64;56.64;56.64];                                          %[m^2];

% Simple orbiting of one body
% options = odeset ('Event',@(t,y) on_ground(t,y,earth_rad));
% [tchas,ychas] = ode45(@(t,y) orbit(t,y,grav_const),[0 7200],[chaser.x0;chaser.xdot0],options);
% [ttar,ytar] = ode45(@(t,y) orbit(t,y,grav_const),[0 7200],[target.x0;target.xdot0],options);


T_on = [randtimes(:,1)*tmax randtimes(:,3)*tmax];
T_off =[randtimes(:,2).*(tmax*ones(3,1)-T_on(:,1)) randtimes(:,4).*(tmax*ones(3,1)-T_on(:,2))] ;

T = [T_on(:,1) T_off(:,1) T_on(:,2) T_off(:,2)];

options = odeset ('Event',@(t,y) on_ground_both(t,y,earth_rad),'MaxStep',10);
[t,y] = ode45(@(t,y) orbit_both(t,y,grav_const,rho,chaser.area(1),target.area(1),chaser.cd(1),target.cd(1),T),[0 tmax],[chaser.x0;target.x0;chaser.xdot0;target.xdot0],options);

% Nice 3D plot

[xsphere,ysphere,zsphere]=sphere(100);
xsphere = earth_rad*xsphere;
ysphere = earth_rad*ysphere;
zsphere = earth_rad*zsphere;

% plot3(ychas(:,1),ychas(:,2),ychas(:,3),'r','LineWidth',3);
% hold on
% plot3(ytar(:,1),ytar(:,2),ytar(:,3),'g','LineWidth',3);
% axis equal
% surf(xsphere,ysphere,zsphere,ones(size(xsphere)))

for i=1:size(t,1)
    hold off
    plot3(y(1:i,1),y(1:i,2),y(1:i,3),'g')%,'LineWidth',3);
    hold on
    plot3(y(1,1),y(1,2),y(1,3),'go')
    plot3(y(end,1),y(end,2),y(end,3),'gx')
    axis equal
    plot3(y(1:i,4),y(1:i,5),y(1:i,6),'r')%,'LineWidth',3);
    plot3(y(1,4),y(1,5),y(1,6),'ro')
    plot3(y(end,4),y(end,5),y(end,6),'rx')
    drawnow
    %surf(xsphere,ysphere,zsphere,ones(size(xsphere)))
    %legend('Chaser','','','Target','','')
end

