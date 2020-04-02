function [f,df]=smooth_cheb_fun(x)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%
% The function to be optimised when using Pascoletti-Serafini scalarisation
% (i.e. minimising alpha)

f = x(1);
df = 0*x;
df(1) = 1;

return