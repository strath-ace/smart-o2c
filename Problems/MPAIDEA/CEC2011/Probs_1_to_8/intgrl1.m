% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function dy = intgrl(t,x,u)

exp(25*pinv(x+2)*x)
exp(25*x*pinv(x+2))
pause
dy = zeros(2,1);    % a column vector
dy(1) = -(2.+u)*(x+0.25)+(x+0.5)*exp(25*pinv(x+2)*x);
dy(2) = 0.5-x-(x+0.5)*exp(25*x/(x+2));
