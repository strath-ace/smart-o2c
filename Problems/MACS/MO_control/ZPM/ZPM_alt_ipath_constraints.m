function val = ZPM_alt_ipath_constraints(x,u,time,static,scales,constants)

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

x = x.*xscale;
u = u.*uscale;                  % NOTE, HERE u ARE qc1,qc2,qc3,qc4!!!
time = time.*tscale;
static = static.*staticscale;

%% eval objective

r = x(1:3);
w = x(4:6);
h = x(7:9);

val = [h'*h-static^2];

%% scaling???

val = val./[constants.hmax^2];

end