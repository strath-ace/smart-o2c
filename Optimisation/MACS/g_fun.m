function g=g_fun(f,lambda,z,zstar)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%
% Tchebychev scalarisation

f2 = (f-z)./(((zstar-z).*((zstar-z)~=0))+((zstar-z)==0));

[~,index] = max(lambda.*abs(f2));
g = lambda(index).*f2(index);

return