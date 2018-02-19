% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [masked, maskedeq] = mask_constraint_max_decomposition(u_to_opt,par_objfun)
% unscales d, u, evaluates the corresponding objfun and multiplies by sign

% d_true = par_objfun.lb_d' + par_objfun.d.*(par_objfun.ub_d' - par_objfun.lb_d');
% DA VERIFICARE 


d_true = par_objfun.d_belief;                                 % for the moment I don't change d

obj = par_objfun.objective;
%map_u_info = par_objfun.map_u_info{obj};

u_aux = par_objfun.u_belief;

% u_aux(par_objfun.vars_to_opt)=u_to_opt;
% u_true = map_affine(u_aux,map_u_info);

%u_true = map_affine(u_to_opt,map_u_info);
u_aux(par_objfun.vars_to_opt)=u_to_opt;
u_true = u_aux;


func = par_objfun.constraint{1}; %par_objfun.objfun{obj};   HARD CODED
par_func = par_objfun.problem_par_objfun{1};

masked = -func(d_true,u_true,par_func);    % -par_objfun.sign*   VERIFICARE 

maskedeq = [];
end