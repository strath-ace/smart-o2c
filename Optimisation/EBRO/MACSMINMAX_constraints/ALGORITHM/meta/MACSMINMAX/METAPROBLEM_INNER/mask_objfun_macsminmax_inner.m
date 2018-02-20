% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [masked] = mask_objfun_macsminmax_inner(u,par_objfun)
% unscales d, u, evaluates the corresponding objfun and multiplies by sign

d_true = par_objfun.lb_d' + par_objfun.d.*(par_objfun.ub_d' - par_objfun.lb_d');

obj = par_objfun.objective;
map_u_info = par_objfun.map_u_info{obj};

u_true = map_affine(u,map_u_info);


func = par_objfun.objfun{obj};
par_func = par_objfun.problem_par_objfun{obj};

masked = -par_objfun.sign*func(d_true,u_true,par_func);

return