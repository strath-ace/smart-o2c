function [decomposition] = max_min_func_decomposition(i, ii, in, position_FE, decomposition, algo_inner, problem_inner)


bpa = 1;
for k = 1:problem_inner.dim
    bpa = bpa*problem_inner.bpa{k,1}(position_FE(k));  % product of the bpa of the coupled components
end

decomposition{1,i-in.num_functions}.FocalElement{1,ii}.bpa = bpa;
decomposition{1,i-in.num_functions}.FocalElement{1,ii}.n_FE = ii;
decomposition{1,i-in.num_functions}.FocalElement{1,ii}.lb = problem_inner.lb;
decomposition{1,i-in.num_functions}.FocalElement{1,ii}.ub = problem_inner.ub;
decomposition{1,i-in.num_functions}.FocalElement{1,ii}.position = position_FE;

%--------------------------------------------------------------------------
% MAX
%--------------------------------------------------------------------------
if in.output == 0 || in.output == 2
    
    
    
    problem_inner.objfun = @mask_objfun_max_decomposition;
    
    % objective and constraints are defined in different functions
    % Function to optimise
    fitnessfcn.obj       = problem_inner.objfun;
    % Function of constraints
    if isempty(problem_inner.par_objfun.constraint{1})
        fitnessfcn.constr    = [];
    else
        fitnessfcn.constr    = @mask_constraint_max_decomposition;
    end
    
    
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
    
    problem_inner.fitnessfcn = fitnessfcn;
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    problem_inner.par_objfun.d = problem_inner.par_objfun.d_belief;
    problem_inner.par_objfun.u = problem_inner.par_objfun.u_belief;
    problem_inner.par_objfun.lb_d = problem_inner.lb;
    problem_inner.par_objfun.ub_d = problem_inner.ub;
    
    
    
    
    [ u_max_to_opt, fmax_to_opt , ~ , ~ ] = algo_inner.optimise(problem_inner,algo_inner.par);
    
    
    decomposition{1,i-in.num_functions}.FocalElement{1,ii}.upper_f = -fmax_to_opt;
    decomposition{1,i-in.num_functions}.FocalElement{1,ii}.upper_u = u_max_to_opt;
    
end


%--------------------------------------------------------------------------
% MAX
%--------------------------------------------------------------------------
if in.output == 1 || in.output == 2
    
    problem_inner.objfun = @mask_objfun_min_decomposition;
    
    
    % objective and constraints are defined in different functions
    % Function to optimise
    fitnessfcn.obj       = problem_inner.objfun;
    % Function of constraints
    if isempty(problem_inner.par_objfun.constraint{1})
        fitnessfcn.constr    = [];
    else
        fitnessfcn.constr    = @mask_constraint_max_decomposition;
    end
    
    
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
    
    problem_inner.fitnessfcn = fitnessfcn;
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    problem_inner.par_objfun.d = problem_inner.par_objfun.d_plausibility;
    problem_inner.par_objfun.u = problem_inner.par_objfun.u_plausibility;
    problem_inner.par_objfun.lb_d = problem_inner.lb;
    problem_inner.par_objfun.ub_d = problem_inner.ub;
    
    
    
    [ u_min_to_opt, fmin_to_opt , ~ , ~ ] = algo_inner.optimise(problem_inner,algo_inner.par);
    
    decomposition{1,i-in.num_functions}.FocalElement{1,ii}.downer_u = u_min_to_opt;
    decomposition{1,i-in.num_functions}.FocalElement{1,ii}.downer_f = fmin_to_opt;
    
    
end

if  in.output ~= 0 && in.output ~= 1 && in.output ~= 2
    print('error')
    
end



global num_maximization_decomposition
num_maximization_decomposition = num_maximization_decomposition +1;

end