function [ x, fval, exitflag, output ] = optimise_idea2_wrapper(problem,par)

if (isfield(par,'initial_population') && ~isempty(par.initial_population))
    initial_population = par.initial_population;
else
    initial_population = lhsu(problem.lb,problem.ub,par.options(4));
end

[memories,nfeval]=idea2(problem.objfun, @c_unconstrained, initial_population, problem.lb, problem.ub, par.options, problem.par_objfun);

x = memories(1,1:problem.dim);
fval = memories(1,problem.dim+1);
exitflag = 1;
output.nfeval = nfeval;

return