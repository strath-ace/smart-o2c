function val = New4ESA_ei_constraints(x0,u0,t0,xf,uf,tf,x,u,time,static,scales,constants)

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

m = x0(7);
S = static(1);

val = [m/S-700 0;];           % wing loading less than 800kg/m^2

%% normalise constraint

val = val./[700];%xscale(5);xscale(5)];

end