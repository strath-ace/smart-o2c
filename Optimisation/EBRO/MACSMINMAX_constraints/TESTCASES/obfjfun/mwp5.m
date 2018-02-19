% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [f] = mwp5(d,u,par)

f = -(d(1)-1)*u(1)-(d(2)-2)*u(2)-(d(3)-1)*u(3)+2*d(1)^2+3*d(2)^2+d(3)^2-u(1)^2-u(2)^2-u(3)^2;
global nfevalglobal
nfevalglobal = nfevalglobal + 1;

end