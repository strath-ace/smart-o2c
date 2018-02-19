% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [out] = TC_1_constraints(d, u, par)


if u(1) + d(1) < 0
    out = abs(u + sum(d));
else
    out = 0;
end

out = sum(out);
return

