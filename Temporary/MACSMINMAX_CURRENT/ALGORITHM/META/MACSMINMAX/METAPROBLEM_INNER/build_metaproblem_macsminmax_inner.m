function [metaproblem] = build_metaproblem_macsminmax_inner(problem)

% for objective function mask
metaproblem.par_objfun.dim_d = problem.dim_d;
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
dim_u =  problem.dim_u;
metaproblem.dim = dim_u;
metaproblem.lb = zeros(1,dim_u);
metaproblem.ub = ones(1,dim_u);

% function
metaproblem.par_objfun.sign = problem.sign_inner;
metaproblem.objfun = @mask_objfun_macsminmax_inner; % depends on u and par_objfun. Needs to specify par_objfun.d and par_objfun.objective

return