function val = New4ESA_ipath_constraints_no_thrust(x,u,time,static,scales,constants)

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
ahat = alpha*180/pi;
S = static(1);
m = static(2);
w_e = constants.omega_e;

%% Compute density, aero forces, gravity

r = constants.Re+h;
g = constants.mu/r^2;
[~,rho,c] = atmo_ISA_smooth(h);
Mach = (v-constants.omega_e*r)/c;

Cl = constants.cl_fun(ahat,Mach);
Cd = constants.cd_fun(ahat,Mach);
L = 0.5*rho*(v-constants.omega_e*r)^2*S*Cl;
D = 0.5*rho*(v-constants.omega_e*r)^2*S*Cd;

vdot = (-D)/m - g*sin(gamma);
gammadot = (+L)/(m*v)*cos(beta) + cos(gamma)*(v/r-g/v);
psidot = (L/(m*v*cos(gamma))*sin(beta) + v/(r*cos(theta))*cos(gamma)*sin(psi)*sin(theta));
        
%% eval constraints

v_imp = (v-constants.omega_e*r)*3.28084;      % convert velocity to imperial units
rho_imp = rho/515.379; % convert density to imperial units

qr = 17700*rho_imp^0.5*(0.0001*v_imp)^3.07;
qa = constants.c0+constants.c1*ahat+constants.c2*ahat^2+constants.c3*ahat^3;
q = qa*qr*11356.538527;     % convert heat flux density back into metric units

totacc2= vdot^2+v^2*(gammadot^2+psidot^2);
  
val = [-h;                                 % h>0km
      -Mach+0.3;
       q-constants.qmax;                        % heat flux below threshol
       0.5*rho*(v-constants.omega_e*r)^2-60e3;                        % dynamic pressure
       totacc2-constants.maxacc^2;               % total acceleration vector below threshold
       -((ahat-35)/15)^8-((Mach-30)/27)^8+1;               % Exclude (alpha,Mach) in the smoothed rectangle where alpha>20deg and Mach<3
     ];

 % scaling

val = val./[xscale(1);1;constants.qmax;60e3;constants.maxacc^2;1000]; 

end