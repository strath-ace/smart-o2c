function [f,df] = feas_only2(x,x0)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%
% General function used when solving pure feasibility problems, returns the
% squared distance between supplied initial guess and current solution

% ensure inputs are column vectors
x = x(:);
x0 = x0(:);

f = (x-x0)'*(x-x0);
df = 2*(x-x0);

end