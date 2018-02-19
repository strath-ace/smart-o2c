% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [max] = evaluate_max(problem, d, algo_inner)

par_objfun.problem = problem;
par_objfun.objective = problem.n_obj;

n_obj_dec = 1;

[meta_algo_inner_max_u] = problem_max_u(problem, n_obj_dec, par_objfun, d, algo_inner);
meta_algo_inner_max_u.par_objfun.flag = 1; %  max

% CONSTRAINTS
fitnessfcn.obj     = meta_algo_inner_max_u.objfun;
if isempty(meta_algo_inner_max_u.constraints)==1
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


    meta_algo_inner_max_u.fitnessfcn = fitnessfcn;
    
% % Flag to 0: objective and constraints are NOT in the same function
% meta_algo_inner_max_u.fitnessfcn.obj_constr = 0;
% % Weights for penalty
% meta_algo_inner_max_u.fitnessfcn.w_ceq = 100;
% meta_algo_inner_max_u.fitnessfcn.w_c   = 100;
%%%%%%%%%%%%%%%%%%%%%%

[ max_u, fmax_u , ~ , ~ ] = meta_algo_inner_max_u.optimise(meta_algo_inner_max_u, meta_algo_inner_max_u.par);

u_true = map_affine(max_u, meta_algo_inner_max_u.par_objfun.map_u_info{meta_algo_inner_max_u.par_objfun.objective});

%%%%%%%%%%%%%%%%%%%%
max.u{1} = u_true;
max.f = -fmax_u;
max.d = d;

end