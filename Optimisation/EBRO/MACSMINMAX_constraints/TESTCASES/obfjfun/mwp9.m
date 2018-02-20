% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [f] = mwp9(d,u,par)

f = min([3-0.2*d(1)+0.3*u(1), 3+0.2*d(1)-0.1*u(1)]);
global nfevalglobal
nfevalglobal = nfevalglobal + 1;

end