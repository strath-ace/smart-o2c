% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [meta_algo_inner_max_u] = problem_max_u(problem, n_obj, par_objfun, d, algo_inner)

[algo_minmin] = map_u_info_minmin(problem, n_obj, problem);
    
    
    meta_algo_inner_max_u = build_metaproblem_max_u(par_objfun.problem, algo_inner);
    
    meta_algo_inner_max_u.par_objfun.map_u_info{par_objfun.objective} = algo_minmin.par_minmax.surrogate{n_obj}.map_info;
    meta_algo_inner_max_u.par_objfun.d = d;
    meta_algo_inner_max_u.par_objfun.objective = n_obj;
    
end