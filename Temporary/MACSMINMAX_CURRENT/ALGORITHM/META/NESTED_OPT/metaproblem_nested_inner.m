function [metaproblem] = metaproblem_nested_inner(problem,algorithm,objective)

% objectives
metaproblem.n_obj = 1;
% minmax_problem.objfun = [@mv1,@mv3];
minmax_problem.objfun_par = problem.objfun_par(objective);

% design variables
minmax_problem.dim_x = 2;
minmax_problem.lb_x = repmat(1,[dim_d,1]);
minmax_problem.ub_x = repmat(5,[dim_d,1]);

% % uncertain variables
% minmax_problem.dim_u = 2;
% minmax_problem.lb_u = {repmat([-5,-3,1],[dim_u,1]),repmat([-5,-3,1],[dim_u,1])};
% minmax_problem.ub_u = {repmat([-4,0,3],[dim_u,1]),repmat([-4,0,3],[dim_u,1])};

