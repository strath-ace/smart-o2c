function dx = Goddard3_phase2_state_equations(x,u,time,static,scales,constants)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%

% states:   x(1) = h
%           x(2) = v
%           x(3) = m
%
% controls: u(1) = T

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

Ts = constants.sigma*x(2)^2*exp(-x(1)/constants.h0)+x(3)*constants.g+x(3)*constants.g/(1+4*constants.c/x(2)+2*constants.c^2/(x(2)^2))*(constants.c^2/(constants.h0*constants.g)*(1+x(2)/constants.c-1-2*constants.c/x(2)));

dx = [x(2);
     (Ts-constants.sigma*x(2)^2*exp(-x(1)/constants.h0))/x(3)-constants.g;       
     -Ts/constants.c];

%% renormalise output

dx = dx*tscale./(xscale);

end