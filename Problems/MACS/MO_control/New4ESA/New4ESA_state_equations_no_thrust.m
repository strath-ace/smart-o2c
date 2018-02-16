function dx = New4ESA_state_equations_no_thrust(x,u,time,static,scales,constants)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------

% initialisation
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
u = u.*uscale;                  % NOTE, HERE u ARE qc1,qc2,qc3,qc4!!!
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
S = static(1);
m = static(2);
w_e = constants.omega_e;

%% Compute aero forces and gravitational acceleration

r = constants.Re+h;
g = constants.mu/r^2;

[~,rho,c] = atmo_ISA_smooth(h);
ahat = 180*alpha/pi;
Mach = (v-constants.omega_e*r)/c;

Cl = constants.cl_fun(ahat,Mach);
Cd = constants.cd_fun(ahat,Mach);
L = 0.5*rho*(v-constants.omega_e*r)^2*S*Cl;
D = 0.5*rho*(v-constants.omega_e*r)^2*S*Cd;

%% the actual dynamics

dx =   [v*sin(gamma);
        v/r*cos(gamma)*sin(psi)/cos(theta);
        v/r*cos(gamma)*cos(psi);
        (-D)/m - g*sin(gamma);
        (+L)/(m*v)*cos(beta) + cos(gamma)*(v/r-g/v);
        (L/(m*v*cos(gamma))*sin(beta) + v/(r*cos(theta))*cos(gamma)*sin(psi)*sin(theta));
        ];

%% renormalise output

dx = dx*tscale./(xscale);

end