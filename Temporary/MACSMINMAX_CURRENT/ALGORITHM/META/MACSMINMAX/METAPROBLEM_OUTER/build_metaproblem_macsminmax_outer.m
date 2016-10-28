function [metaproblem] = build_metaproblem_macsminmax_outer(problem_fix_d, local_search)

% chromosome
dim_d =  problem_fix_d.par_objfun.dim_d;
metaproblem.dim = dim_d;
metaproblem.lb = zeros(1,dim_d);
metaproblem.ub = ones(1,dim_d);

% function
metaproblem.par_objfun.problem_fix_d = problem_fix_d;
metaproblem.par_objfun.u_record = [];
metaproblem.par_objfun.local_search = local_search;
metaproblem.par_objfun.objectives = 1:problem_fix_d.n_obj;
metaproblem.objfun = @mask_objfun_macsminmax_outer; % depends on d and par_objfun. Needs to specify par_objfun.(all that is above)

return