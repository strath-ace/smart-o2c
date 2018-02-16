function val = betts_descent_ipath_constraints(x,u,time,static,scales,constants)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%

% states:   x(1) = h (altitude) [ft]
%           x(2) = phi (longitude) [rad]
%           x(3) = theta (latitude) [rad]
%           x(4) = v (velocity) [ft/s]
%           x(5) = gamma (fpa) [rad]
%           x(6) = psi (azimuth) [rad]
%
% controls: u(1) = alpha (aoa) [rad]
%           u(2) = beta (bank) [rad]

%% constants

xscale = scales.xscale;
uscale = scales.uscale;
tscale = scales.tscale;
staticscale = scales.other_vars_scale;

%% de normalising input

x = x.*xscale;
u = u.*uscale;                 
time = time.*tscale;
static = static.*staticscale;

%% rename variables

h = x(1);
theta = x(3);
v = x(4);
gamma = x(5);
psi = x(6);
qu = static(1);
alpha = u(1);
beta = u(2);
ahat = alpha*180/pi;
r = constants.Re+h;

%% Compute density, aero forces, gravity

rho = constants.rho0*exp(-h/constants.hr);
Cl = constants.a0+constants.a1*ahat;
L = 0.5*rho*v^2*constants.S*Cl;
g = constants.mu/r^2;

%% eval constraints

qr = 17700*rho^0.5*(0.0001*v)^3.07;
qa = constants.c0+constants.c1*ahat+constants.c2*ahat^2+constants.c3*ahat^3;
gammadot = L/(constants.m*v)*cos(beta) + cos(gamma)*(v/r-g/v);
psidot = L/(constants.m*v*cos(gamma))*sin(beta) + v/(r*cos(theta))*cos(gamma)*sin(psi)*sin(theta);

val = [(qr*qa)-qu;
       gammadot-constants.zeta;
       psidot-constants.zeta;
       -gammadot-constants.zeta;
       -psidot-constants.zeta;
       ];

%% scaling???

val = val./[1000;0.2;0.2;0.2;0.2]; % 1000 is maximum for parameter qu

end