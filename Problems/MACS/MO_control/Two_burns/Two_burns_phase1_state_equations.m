function dx = Two_burns_phase1_state_equations(x,u,time,static,scales,constants)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%

% states:   x(1) = p
%           x(2) = f
%           x(3) = g
%           x(4) = h
%           x(5) = k
%           x(6) = L

%
% controls: u(1) = uy

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

%% the actual dynamics

q = 1+x(2)*cos(x(6))+x(3)*sin(x(6));

dx = [0;
      0;
      0;
      0;
      0;
      (constants.mu*x(1))^0.5*(q/x(1))^2];

%% renormalise output

dx = dx*tscale./(xscale);

end