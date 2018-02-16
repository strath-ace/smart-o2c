function dx =Docking_state_equations(x,u,time,static,scales,constants)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%

% states: x, y, z, vx, vy, vz, q1, q2, q3, q4, w1, w2, w3, p1, p2, p3, p4, phi1, phi2, phi3
% controls: alphax, alphay, alphaz, taux, tauy, tauz

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

%% Preliminaries

q = x(7:10);
q = q/(q'*q)^0.5;  % forceful normalisation to unit quaternion
omega = x(11:13);
p = x(14:17);
p = p/(p'*p)^0.5;  % forceful normalisation to unit quaternion
phi = x(18:20);

Omega = Omega_fun(omega(1),omega(2),omega(3));
Phi = Omega_fun(phi(1),phi(2),phi(3));

%% the actual dynamics

dx = [x(4);
      x(5);
      x(6);
      2*constants.n*x(5)+3*constants.n^2*x(1)+u(1)/constants.m;
      -2*constants.n*x(4)+u(2)/constants.m;
      -constants.n^2*x(3)+u(3)/constants.m;
      0.5*Omega*q;
      (omega(2)*omega(3)*(constants.J22-constants.J33)+u(4))/constants.J11;
      (omega(1)*omega(3)*(constants.J33-constants.J11)+u(5))/constants.J22;
      (omega(1)*omega(2)*(constants.J11-constants.J22)+u(6))/constants.J33;
      0.5*Phi*p;
      (phi(2)*phi(3)*(constants.K22-constants.K33))/constants.K11;
      (phi(1)*phi(3)*(constants.K33-constants.K11))/constants.K22;
      (phi(1)*phi(2)*(constants.K11-constants.K22))/constants.K33];
  
%% renormalise output

dx = dx*tscale./(xscale);

end