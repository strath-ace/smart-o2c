function [ x, fval, exitflag, output ] = optimise_constraint(problem,par)
% Maximisation of the constraint in the uncertain space.
% 
% min(-C)
if (isfield(par,'initial_population') && ~isempty(par.initial_population))
    initial_population = par.initial_population;
else
    % Initialise populations
    initial_population = zeros(par.n_agents,problem.dim,par.n_populations);

    for s = 1 : par.n_populations
        % pop = lhsdesign(par.n_agents,problem.dim,'criterion','maximin').*repmat(problem.ub-problem.lb,par.n_agents,1)+repmat(problem.lb,par.n_agents,1);
        pop = lhsgen(par.n_agents,problem.dim).*repmat(problem.ub-problem.lb,par.n_agents,1)+repmat(problem.lb,par.n_agents,1);

        
        initial_population(:,:,s) = pop;
    end
end


% Run MP-AIDEA


% %% %%%%%%%%%%%%%%%%%%%
% problem.par_objfun.sign_constraint  = -1;
% %% %%%%%%%%%%%%%%%%%%%

% options = optimset('Display','none','MaxFunEvals',50*100,'TolFun',1e-8,...%'LargeScale','off',...
% 'Algorithm','sqp'); % add a converged stop condition. in the original there was one but wrongly implemented
% 
% [x,fval,~,output] = fmincon(@mask_constraints_macsminmax_max_constraint,initial_population,[],[],[],[],problem.lb, problem.ub, [],options, problem.par_objfun);

% Function to optimise
fitnessfcn.obj       = @mask_constraints_macsminmax_max_constraint;
% Function of constraints
fitnessfcn.constr    = [];%problem.par_objfun.mask_constraints;


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



% Flag to 0: objective and constraints are NOT in the same function
% fitnessfcn.obj_constr = 0;
% Weights for penalty
% fitnessfcn.w_ceq = [];
% fitnessfcn.w_c = [];


[memories_record, memories, archivebest, population_evolution, vval_evolution,...
    B_mean, delta_local, inite, iglob, options, exitflag] = MP_AIDEA(fitnessfcn, problem.lb, problem.ub, initial_population, par, problem.par_objfun);

%%%%%%%%%%%%%%%%%%%%%
% problem.par_objfun.sign_constraint  = 1;
%%%%%%%%%%%%%%%%%%%%%

% Output: minima and minima's objective value
x    = zeros(size(memories_record,3), problem.dim);
fval = zeros(size(memories_record,3), 1);

for i = 1 : size(memories_record,3)
    x(i,:)  = memories_record(end, 1:end-1, i);
    fval(i) = - memories_record(end, end, i);
end

output.memories             = memories;
output.archivebest          = archivebest;
output.population_evolution = population_evolution;
output.vval_evolution       = vval_evolution;
output.B_mean               = B_mean;
output.delta_local          = delta_local;
output.number_LR            = inite;
output.number_GR            = iglob;
output.nfeval               = options.nFeValMax;

return