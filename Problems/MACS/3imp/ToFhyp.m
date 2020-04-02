% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function t12 = ToFhyp(F1,F2,amax,ecc)

amax = min(amax,0);

t12 = sqrt(-amax^3)*(ecc*(sinh(F2)-sinh(F1)) - F2 + F1);

