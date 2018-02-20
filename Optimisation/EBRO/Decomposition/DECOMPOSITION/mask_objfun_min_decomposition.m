% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [masked] = mask_objfun_min_decomposition(u_to_opt,par_objfun)

d_true = par_objfun.d_plausibility;                                 % for the moment I don't change d

obj = par_objfun.objective;

u_aux = par_objfun.u_plausibility; 

u_aux(par_objfun.vars_to_opt)=u_to_opt;

u_true = u_aux;

func = par_objfun.objfun{obj};

par_func = par_objfun.problem_par_objfun{1};

masked = func(d_true,u_true,par_func);    

end