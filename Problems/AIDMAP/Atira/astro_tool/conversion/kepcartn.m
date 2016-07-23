function [car]=kepcartn(kep, V, gam0);
%
%   [car]=kepcart(kep, V);
%  
%   convertes from keplerian orbital elements to Cartesian coordinates
%   
%   INPUT
%               kep = vector of keplerian elements
%                          a, e, i, Omega, omega, tho
%                          angles are in degrees
%               mu  = planetary gravity constat
%
%   OUTPUT
%                car = cartesian coordinates
%                         the first three components are the position
%                         vector while the last 3 are the velocity vector
%
%  (c) Massimiliano Vasile 2002
%

a   = kep(1);
e   = kep(2);
i   = kep(3);
Om  = kep(4);
om  = kep(5);
tho = kep(6);

i=i*pi/180;
Om=Om*pi/180;
om=om*pi/180;
tho=tho*pi/180;

om = om+tho;

% Rotation matrix
R(1)=cos(om)*cos(Om)-sin(om)*cos(i)*sin(Om);
R(2)=cos(om)*sin(Om)+sin(om)*cos(i)*cos(Om);
R(3)=sin(om)*sin(i);

gam= atan2(e*sin(tho), 1+e*cos(tho))+gam0*pi/180;

% R(1, 2)=-sin(om)*cos(Om)-cos(om)*cos(i)*sin(Om);
% R(2, 2)=-sin(om)*sin(Om)+cos(om)*cos(i)*cos(Om);
% R(3, 2)=cos(om)*sin(i);
% 
% R(1, 3)=sin(i)*sin(Om);
% R(2, 3)=-sin(i)*cos(Om);
% R(3, 3)=cos(i);


%  In plane     Parameters
% p = a*(1-e^2);
% r=p/(1+e*cos(tho));
% xp=r*cos(tho);
% yp=r*sin(tho);
% wom_dot = sqrt(mu*p)/r^2;
% r_dot      = sqrt(mu/p)*e*sin(tho);
% vxp=r_dot*cos(tho)-r*sin(tho)*wom_dot;
% vyp=r_dot*sin(tho)+r*cos(tho)*wom_dot;

% 3D cartesian vector
% car(1)=R(1, 1)*xp+R(1, 2)*yp;
% car(2)=R(2, 1)*xp+R(2, 2)*yp;
% car(3)=R(3, 1)*xp+R(3, 2)*yp;

car(1)=R(1)*sin(gam) - cos(gam)*(cos(Om)*sin(om)+sin(Om)*cos(i)*cos(om));
car(2)=R(2)*sin(gam) - cos(gam)*(sin(Om)*sin(om)-cos(Om)*cos(i)*cos(om));
car(3)=R(3)*sin(gam) + cos(gam)*cos(om)*sin(i);

car(1:3)=car(1:3)*V;
return
