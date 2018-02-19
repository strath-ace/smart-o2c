% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [masked] = mask_objfun_ideaminmax_s_inner(u,par_objfun)

n_FE = par_objfun.surrogate{par_objfun.objective}.find([par_objfun.d u],  par_objfun.surrogate{par_objfun.objective});
[y_aux,mse_aux] = par_objfun.surrogate{par_objfun.objective}.predictor([par_objfun.d u], par_objfun.surrogate{par_objfun.objective}.model{n_FE});
y = -par_objfun.sign*y_aux;
mse = mse_aux;

masked = -par_objfun.surrogate{par_objfun.objective}.indicator_u(y,mse,par_objfun.ymin(par_objfun.objective)); % negative because we will always maximize PI or EI

return