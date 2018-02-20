% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [f] = mv3(d,u,par)
% Multidimensional test function MV3

f = sum( (5-d).*(1+cos(u)) + (d-1).*(1+sin(u)) );

global nfevalglobal
nfevalglobal = nfevalglobal + 1;

end
