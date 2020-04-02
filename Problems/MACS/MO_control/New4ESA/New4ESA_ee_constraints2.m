function val = New4ESA_ee_constraints2(x0,u0,t0,xf,uf,tf,x,u,time,static,scales,constants)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------


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

S = static(1);
m0 = static(2);

% wing mass estimation
Aexp = 0.7283*S;
Sbody = 0.6970*S;
bstr = 1.486.*S.^0.5;
bbody = 0.2785.*S.^0.5;
troot = 0.087.*S.^0.5;

Nz = 1.5*3; % load factor 3, factor of safety 1.5
eta = 0.15; % control configured vehicle
Kwing= 0.214; % Organic composite honeycomb, no TPS (the lightest, also is the only composite and X33 was all composite)
Kct = 0.0267; % dry carry-thru (integral, the lightest(


Mwing = (Nz.*m0*1./(1+eta*(Sbody/Aexp))).^0.386.*(Aexp./troot).^0.572.*(Kwing.*bstr.^0.572+Kct*bbody.^0.572);


val = [ xf(7)-(10e3+Mwing) 0;% final mass equals static mass of 10ton + wings; 
       tf-2*t0 0];           % this should force elemets of phase 1 and 2 to be equally spaced

%% normalise constraint

val = val./[xscale(7);tscale];

end