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

% CONSTRAINTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
metaproblem.par_objfun.constraint{1} = problem.constraints{1};
if isempty(metaproblem.par_objfun.constraint{1})
    metaproblem.par_objfun.mask_constraints = [];
else
    metaproblem.par_objfun.mask_constraints = @mask_constraints_macsminmax_inner;  
end
% metaproblem.par_objfun.constraints_flag = 1; % do constraints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% chromosome
dim_u =  problem.dim_u;
metaproblem.dim = dim_u;
metaproblem.lb = zeros(1,dim_u);
metaproblem.ub = ones(1,dim_u);

% function
metaproblem.par_objfun.sign = problem.sign_inner;
metaproblem.objfun = @mask_objfun_macsminmax_inner; % depends on u and par_objfun. Needs to specify par_objfun.d and par_objfun.objective


% objective and constraints are defined in different functions

% Function to optimise
fitnessfcn.obj       = metaproblem.objfun;
% Function of constraints
fitnessfcn.constr    = metaproblem.par_objfun.mask_constraints;
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

metaproblem.fitnessfcn = fitnessfcn;


% Flag to 0: objective and constraints are NOT in the same function
% metaproblem.fitnessfcn.obj_constr = 0;
% Weights for penalty
% metaproblem.fitnessfcn.w_ceq = 100;
% metaproblem.fitnessfcn.w_c   = 100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return