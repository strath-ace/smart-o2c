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

% CONSTRAINTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
metaproblem.par_objfun.constraint{1} = problem_fix_d.par_objfun.constraint;
if isempty(problem_fix_d.par_objfun.constraint{1})
    metaproblem.par_objfun.mask_constraints = [];
else
    metaproblem.par_objfun.mask_constraints = @mask_constraints_macsminmax_outer; %[]; %
end
% metaproblem.par_objfun.constraints_flag = 0;


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
% metaproblem.fitnessfcn.w_c = 100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


return