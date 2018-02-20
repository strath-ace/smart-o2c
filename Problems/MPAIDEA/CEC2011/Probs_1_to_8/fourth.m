% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function f=fourth(x,u) 
[ps d]=size(x)
f=sum(Y.^2,2)+0.1*repmat(u*u,ps,1);
% f(i)=x(1)^2+x(2)^2+0.1*u^2;
