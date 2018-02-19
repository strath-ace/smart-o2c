% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [metaproblem] = build_metaproblem_ideaminmax_s_outer(problem)

% chromosome
dim_d =  problem.dim_d;
metaproblem.dim = dim_d;
metaproblem.lb = zeros(1,dim_d);
metaproblem.ub = ones(1,dim_d);

metaproblem.par_objfun.objectives = 1:problem.n_obj;
metaproblem.par_objfun.surrogate = [];
metaproblem.par_objfun.indicator_d = [];
metaproblem.par_objfun.ymin = [];
metaproblem.par_objfun.sign_inner = problem.sign_inner;
% function
metaproblem.objfun = @mask_objfun_ideaminmax_s_outer; % depends on d and par_objfun. Needs to specify par_objfun.surrogate.model, par_objfun.ymin
return