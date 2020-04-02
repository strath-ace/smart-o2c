function val = ZPM_alt_eepath_constraints(x0,u0,t0,xf,uf,tf,x,u,time,static,scales,constants)

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

r = xf(1:3);
w = xf(4:6);
h = xf(7:9);

C = eye(3)+2/(1+r'*r)*(constants.vcross(r)*constants.vcross(r)-constants.vcross(r));   % Rotation matrix
C2 = C(:,2);
C3 = C(:,3);

w0 = -constants.n*C2;
taud = 3*(constants.n)^2*constants.vcross(C3)*constants.J*C3;    % only gravity gradient torque, no aero for now

rdot = 0.5*(r*r'+eye(3)+constants.vcross(r))*(w-w0);
omegadot = constants.Jinv*(taud - constants.vcross(w)*(constants.J*w+h));

val = [rdot zeros(3,1); omegadot zeros(3,1)]; 

end