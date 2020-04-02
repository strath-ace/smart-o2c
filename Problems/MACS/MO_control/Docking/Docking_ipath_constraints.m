function val = Docking_ipath_constraints(x,u,time,static,scales,constants)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%

% states:   x(1) = x
%           x(2) = l
%           x(3) = z
%
% controls: u(1) = u

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

%% the actual constraint

%q = x(7:10);
%Q = Q_fun(q(1),q(2),q(3),q(4));

val = [% u(1:3)'*u(1:3)-constants.umax;
       -(x(1:3)'*x(1:3))^0.5+2];

%% renormalise output

val = val./[(xscale(1:3)'*xscale(1:3))^0.5];
%val = val./[constants.umax;(xscale(1:3)'*xscale(1:3))^0.5];

end