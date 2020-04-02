function val = ZPM_alt_objective(x0,u0,t0,xf,uf,tf,x,u,time,static,scales,constants)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%

% states:   x(1) = r1
%           x(2) = r2
%           x(3) = r3
%           x(4) = w1
%           x(5) = w2
%           x(6) = w3
%           x(7) = h1
%           x(8) = h2
%           x(9) = h3
%
% controls: u(1) = ux
%           u(2) = uy
%           u(3) = uz

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

%% eval objective

r = x(1:3);
w = x(4:6);
h = x(7:9);

val = [static(1) 0; 0 1e-6*(u'*u)];

end