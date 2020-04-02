function val = New4ESA_ee_constraints(x0,u0,t0,xf,uf,tf,x,u,time,static,scales,constants)

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

x0 = x0.*xscale;
u0 = u0.*uscale;
t0 = t0.*tscale;
xf = xf.*xscale;
uf = uf.*uscale;
tf = tf.*tscale;
x = x.*xscale;
u = u.*uscale;                  
time = time.*tscale;
static = static.*staticscale;

%% eval constraint

h = xf(1);
phi = xf(2);
theta = xf(3);
v = xf(4);
gamma = xf(5);
psi = xf(6);
alpha = uf(1);
beta = uf(2);
S = static(1);
m = static(2);
w_e = constants.omega_e;

r = constants.Re+h;
g = constants.mu/r^2;

[~,rho,c] = atmo_ISA_smooth(h);
ahat = 180*alpha/pi;
Mach = (v-constants.omega_e*r)/c;

Cl = constants.cl_fun(ahat,Mach);
Cd = constants.cd_fun(ahat,Mach);
L = 0.5*rho*v^2*S*Cl;
D = 0.5*rho*v^2*S*Cd;

val = [% L/(m*v*cos(gamma))*sin(beta) - v/r*cos(gamma)*cos(psi)*tan(theta)+2*w_e*(sin(psi)*cos(theta)*tan(gamma)-sin(theta))-w_e^2*r/(v*cos(gamma))*cos(theta)*sin(gamma)*cos(psi) 0
       Mach-0.4 0];  

%% normalise constraint

val = val./[0.4];

end