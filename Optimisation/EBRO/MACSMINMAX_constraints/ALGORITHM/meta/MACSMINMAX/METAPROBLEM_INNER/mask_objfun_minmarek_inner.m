% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [masked] = mask_objfun_minmarek_inner(u,par_objfun)

[y,mse] = par_objfun.surrogate.predictor(u, par_objfun.surrogate.model);
masked = -par_objfun.surrogate.indicator(y,mse,par_objfun.ymin); % negative because we will always maximize PI or EI

return