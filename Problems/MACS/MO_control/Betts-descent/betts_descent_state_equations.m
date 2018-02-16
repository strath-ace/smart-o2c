function dx = betts_descent_state_equations(x,u,time,static,scales,constants)

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
phi = x(2);
theta = x(3);
v = x(4);
gamma = x(5);
psi = x(6);

alpha = u(1);
beta = u(2);

%% Compute aero forces and gravitational acceleration

r = constants.Re+h;
g = constants.mu/r^2;
rho = constants.rho0*exp(-h/constants.hr);
ahat = 180*alpha/pi;
Cl = constants.a0+constants.a1*ahat;
Cd = constants.b0+constants.b1*ahat+constants.b2*ahat^2;
L = 0.5*rho*v^2*constants.S*Cl;
D = 0.5*rho*v^2*constants.S*Cd;

%% the actual dynamics

dx =   [v*sin(gamma);
        v/r*cos(gamma)*sin(psi)/cos(theta);
        v/r*cos(gamma)*cos(psi);
        -D/constants.m - g*sin(gamma);
        L/(constants.m*v)*cos(beta) + cos(gamma)*(v/r-g/v);
        L/(constants.m*v*cos(gamma))*sin(beta) + v/(r*cos(theta))*cos(gamma)*sin(psi)*sin(theta);
        ];

%% renormalise output

dx = dx*tscale./(xscale);

end