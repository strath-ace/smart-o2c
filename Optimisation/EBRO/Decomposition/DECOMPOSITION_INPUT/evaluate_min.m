% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [min] = evaluate_min(problem, d, algo_inner)

par_objfun.problem = problem;
par_objfun.objective = problem.n_obj;


n_obj_dec = 1;

par_objfun.problem.flag = -1;  % min
par_objfun.problem.sign = -1;


% CONSTRAINTS
fitnessfcn.obj     = @meta_mask_objfun_min_u;
if isempty(problem.constraints{1, 1})==1
    fitnessfcn.constr = [];
else    
    fitnessfcn.constr = @metamask_constraints; %meta_algo_inner_max_u.constraints;
end
% Flag to 0: objective and constraints are NOT in the same function
fitnessfcn.obj_constr = 0;
% How to handle constraints: set to 1 for weighted constraints with fixed
% weights, or to 0 for penalty with no weights
fitnessfcn.weighted = 0;
% If the constraints are handled without weights, then define a tolerance
% for the violation of the equality constraints
fitnessfcn.ceq_eps = 1e-6;
% Weights for penalty if fitnessfcn.weighted == 1
fitnessfcn.w_ceq = 100;
fitnessfcn.w_c = 100;


    



[meta_algo_inner_max_u] = problem_max_u(problem, n_obj_dec, par_objfun, d, algo_inner);
meta_algo_inner_max_u.par_objfun.flag = -1; %  min
meta_algo_inner_max_u.fitnessfcn = fitnessfcn;


[ min_u, fmin_u , ~ , ~ ] = meta_algo_inner_max_u.optimise(meta_algo_inner_max_u, meta_algo_inner_max_u.par);

u_true_min = map_affine(min_u, meta_algo_inner_max_u.par_objfun.map_u_info{meta_algo_inner_max_u.par_objfun.objective});




%%%%%%%%%%%%%%%%%
min.f = fmin_u;
min.u{1} = u_true_min;
min.d = d;

end