function dx = Free_robot_state_equations(x,u,time,static,scales,constants)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%

% states:   x(1) = x
%           x(2) = vx
%           x(3) = y
%           x(4) = vy
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

%% the actual dynamics

dx = [x(4);
      x(5);
      x(6);
      (u(1)-u(2)+u(3)-u(4))*cos(x(3));
      (u(1)-u(2)+u(3)-u(4))*sin(x(3));
      constants.alpha*(u(1)-u(2))-constants.alpha*(u(3)-u(4))];

%% renormalise output

dx = dx*tscale./(xscale);

end