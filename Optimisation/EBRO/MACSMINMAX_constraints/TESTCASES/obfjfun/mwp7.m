% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [f] = mwp7(d,u,par)

f = 2*d(1)*d(5)+3*d(4)*d(2)+d(5)*d(3)+5*d(4)^2+5*d(5)^2-d(4)*(u(4)-u(5)-5)+d(5)*(u(4)-u(5)+3)+sum(u(1:3).*(d(1:3).^2-1))-sum(u.^2);
global nfevalglobal
nfevalglobal = nfevalglobal + 1;

end