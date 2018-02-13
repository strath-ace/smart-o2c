function [metaproblem] = build_metaproblem_so_minmin(problem)

% for objective function mask
dim_d = problem.dim_d;
metaproblem.par_objfun.dim_d = dim_d;

metaproblem.par_objfun.lb_d = problem.lb_d;
metaproblem.par_objfun.ub_d = problem.ub_d;
metaproblem.n_obj = problem.n_obj;

for obj = 1:problem.n_obj
    metaproblem.par_objfun.objfun{obj} = problem.objfun{obj};
    metaproblem.par_objfun.problem_par_objfun{obj} = problem.par_objfun{obj};
    % metaproblem.par_objfun.lb_u{obj} = problem.lb_u{obj};
    % metaproblem.par_objfun.ub_u{obj} = problem.ub_u{obj};
    metaproblem.par_objfun.map_u_info{obj} = get_map_info(problem.lb_u{obj}, problem.ub_u{obj});
end

% chromosome
dim_u = problem.dim_u;
metaproblem.par_objfun.dim_u = dim_u;

metaproblem.dim = dim_u+dim_d;
metaproblem.lb = zeros(1,metaproblem.dim);
metaproblem.ub = ones(1,metaproblem.dim);

% function
metaproblem.objfun = @mask_objfun_so_minmin; % depends on u and par_objfun. Needs to specify par_objfun.d and par_objfun.objective





% CONSTRAINTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
metaproblem.par_objfun.constraint{1} = problem.constraints{1};
metaproblem.par_objfun.mask_constraints = @mask_constraints_macsminmax_inner; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




return