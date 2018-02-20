% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [Cud,Cd,Cu,CEQud,CEQd,CEQu] = c_unconstrained(x_u,x_d,l,func,vub_d,vlb_d,bin,bout,id_par,n,par,flag)

% Inequality and equality constraints. Constraint functions can be
% nonlinear. The constraints are in the form
%   C(x) <= 0
%   CEQ(x) = 0
%
% Example:
%   x(1)^2 + x(2) <= 0.5
%   x(1)*x(2) = 0.2
% becomes
%   C = x(1)^2 + x(2) - 0.5;
%   CEQ = x(1)*x(2) - 0.2;
%
% The suffix refers to:
%   ud: constraints involving both u and d
%    d: constraints involving d only
%    u: constraints involving u only
%
% Added heading, Simone Alicino, 5/11/2012
% Corrected input and added conversion, Simone Alicino, 8/11/2012

% x_d = x_d(:);
% x_d = x_d';
% x(id_par.d) = x_d;
% x(id_par.u) = x_u;

% Conversion from unit hypercube to actual space.
% d = vlb_d + x_d.*(vub_d - vlb_d);
% u = affine_int(x_u,bout,bin);
% x(id_par.d) = d;
% x(id_par.u) = u;

% Constraints initialization
Cud = 0;
Cd = 0;
Cu = 0;
CEQud = 0;
CEQd = 0;
CEQu = 0;

% Inequality constraints
% Cu(1) = 0.1 - x(id_par.u(1));    % u1 > 0.1
% Cu(2) = x(id_par.u(2)) - 0.8;    % u2 <= 0.8
% 
% Cd = 0.1 - x(id_par.d(1));       % d1 > 0.1
% 
% Cud = 0.03 - x(id_par.d(1))*x(id_par.u(1));   % d1*u1 > 0.03

% Equality constraints


return
