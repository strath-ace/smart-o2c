function val = link_ph_2_3(phase_a,phase_b)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%

%% constants

xscale_a = phase_a.scales.xscale;
uscale_a = phase_a.scales.uscale;
tscale_a = phase_a.scales.tscale;
staticscale_a = phase_a.scales.other_vars_scale;
constants_a = phase_a.constants;

xscale_b = phase_b.scales.xscale;
uscale_b = phase_b.scales.uscale;
tscale_b = phase_b.scales.tscale;
staticscale_b = phase_b.scales.other_vars_scale;
constants_b = phase_b.constants;

%% de normalising input

% phase a

x0_a = phase_a.x0.*xscale_a;
u0_a = phase_a.u0.*uscale_a;                  
t0_a = phase_a.t0.*tscale_a;
xf_a = phase_a.xf.*xscale_a;
uf_a = phase_a.uf.*uscale_a;                  
tf_a = phase_a.tf.*tscale_a;
static_a = phase_a.static.*staticscale_a;

% phase b

x0_b = phase_b.x0.*xscale_b;
u0_b = phase_b.u0.*uscale_b;                  
t0_b = phase_b.t0.*tscale_b;
xf_b = phase_b.xf.*xscale_b;
uf_b = phase_b.uf.*uscale_b;                  
tf_b = phase_b.tf.*tscale_b;
static_b = phase_b.static.*staticscale_b;

%% Linkage function

val = [xf_a(1:6)-x0_b(1:6);tf_a-t0_b];         

%% normalise output

val = val./[xscale_b(1:6);tscale_b];

end
