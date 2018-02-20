% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [metaproblem] = build_metaproblem_ideaminmax_s_inner(problem)




% chromosome
dim_u =  problem.dim_u;  % length(problem.u_max_tot); %
metaproblem.dim = dim_u;
metaproblem.lb = zeros(1,dim_u);
metaproblem.ub = ones(1,dim_u);

%%%
for obj = 1:problem.n_obj
    metaproblem.par_objfun.objfun{obj} = problem.objfun{obj};
    metaproblem.par_objfun.constraint{obj} = problem.constraints{obj};
    metaproblem.par_objfun.problem_par_objfun{obj} = problem.par_objfun{obj};
    % metaproblem.par_objfun.lb_u{obj} = problem.lb_u{obj};
    % metaproblem.par_objfun.ub_u{obj} = problem.ub_u{obj};
    metaproblem.par_objfun.map_u_info{obj} = get_map_info(problem.lb_u{obj}, problem.ub_u{obj});
end
%%%


metaproblem.par_objfun.surrogate = [];

metaproblem.par_objfun.sign = problem.sign_inner; %dunno if necessary
metaproblem.par_objfun.ymin = [];

metaproblem.objfun = @mask_objfun_ideaminmax_s_inner; %depends on u and par_objfun. Needs to specify par_objfun.surrogate.model, par_objfun.ymin

if isempty(problem.constraints{1})
    metaproblem.mask_constraints = [];
else
    metaproblem.mask_constraints = @mask_constraints_macsminmax_inner; 
end



return