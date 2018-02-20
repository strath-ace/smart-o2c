function val = Docking_bc(x0,u0,t0,xf,uf,tf,x,u,time,static,scales,constants)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%----Copyright (C) 2018 University of Strathclyde and Authors-----------

% states:   x(1) = r
%           x(2) = vr
%           x(3) = theta
%           x(4) = vt
%
% controls: u(1) = u

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

q = xf(7:10);
omega = xf(11:13);
p = xf(14:17);
phi = xf(18:20);
S = Q_fun(q(1),q(2),q(3),q(4))';
T = Q_fun(p(1),p(2),p(3),p(4))';

% CONDITIONS FROM MICHAELS ET AL, MORE CORRECT AS THEY (IMPLICITLY) IMPOSE
% MATCHING OF FINAL RELATIVE ATTITUDE AND ANGULAR VELOCITIES BY IMPOSING 
% S=T AT FINAL TIME AND THUS REWRITING THE CONDITIONS AS:

% val = [xf(1:3)-R*(constants.bvec-constants.avec) zeros(3,1);
%        xf(4:6)-cross(R*omega,R*(constants.bvec-constants.avec)) zeros(3,1);
%        q-p zeros(4,1);
%        omega-phi zeros(3,1)];

% CONDITIONS FROM BETTS, MISSING RELATIVE ATTITUDE AND RELATIVE ANGULAR VELOCITIES!!

val = [xf(1:3)+S*constants.avec-T*constants.bvec zeros(3,1)                                 % position of docking points
      xf(4:6)+cross(S*omega,S*constants.avec)-cross(T*phi,T*constants.bvec) zeros(3,1) ];  % linear velocity of docking points
                                                                                            

%% normalise constraint

val = val./[xscale(1:6) ones(6,1)];
%val = val./[xscale(1:13) ones(13,1)];

end