% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [problem_inner] = build_problem_inner(problem, minmax, minmin, n_obj)


dim_u =  problem.dim_u;
problem_inner.dim = dim_u;
problem_inner.lb = zeros(1,dim_u);
problem_inner.ub = ones(1,dim_u);
problem_inner.par_objfun.sign = problem.sign_inner; 
problem_inner.par_objfun.ymin = [];

problem_inner.par_objfun.objective = n_obj;

if isfield(minmax,'d')
    problem_inner.par_objfun.d_belief = minmax.d;
    problem_inner.par_objfun.u_belief = minmax.u;
end

if isfield(minmin,'d')
    problem_inner.par_objfun.d_plausibility  = minmin.d;
    problem_inner.par_objfun.u_plausibility = minmin.u;  
end

problem_inner.par_objfun.objfun = problem.objfun;          
problem_inner.par_objfun.constraint = problem.constraints;

problem_inner.par_objfun.problem_par_objfun{1} = problem.par_objfun{1};

return