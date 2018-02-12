function [ x, fval, exitflag, output ] = optimise_mpaidea_wrapper(problem,par)

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



% %% objective and constraints are defined in different functions
% 
% % Function to optimise
% fitnessfcn.obj       = problem.objfun;
% % Function of constraints
% fitnessfcn.constr    = problem.par_objfun.mask_constraints;
% % Flag to 0: objective and constraints are NOT in the same function
% fitnessfcn.obj_constr = 0;
% % Weights for penalty
% fitnessfcn.w_ceq = 100;
% fitnessfcn.w_c = 100;

fitnessfcn = problem.fitnessfcn;

% MP-AIDEA optimisation
% [x,fval,exitflag,output] = optimise_mpaidea(fitnessfcn, problem.lb, problem.ub, options);
%%

% Run MP-AIDEA
[memories_record, memories, archivebest, population_evolution, vval_evolution,...
    B_mean, delta_local, inite, iglob, options, exitflag] = MP_AIDEA(fitnessfcn, problem.lb, problem.ub, initial_population, par, problem.par_objfun);

% Output: minima and minima's objective value
x    = zeros(size(memories_record,3), problem.dim);
fval = zeros(size(memories_record,3), 1);

for i = 1 : size(memories_record,3)
    x(i,:)  = memories_record(end, 1:end-1, i);
    fval(i) = memories_record(end, end, i);
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