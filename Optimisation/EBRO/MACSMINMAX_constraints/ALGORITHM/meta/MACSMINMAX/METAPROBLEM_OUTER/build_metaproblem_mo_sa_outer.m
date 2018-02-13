function [metaproblem] = build_metaproblem_mo_sa_outer(problem)

% chromosome
dim_d =  problem.dim_d;
metaproblem.dim = dim_d;
metaproblem.lb = zeros(1,dim_d);
metaproblem.ub = ones(1,dim_d);

% surrogate
% These are hardcoded (HC) at the moment but should really come from par_minmax for flexibility
    surrogate.method = 'kriging';
    surrogate.corrfun = @corrgauss;
    surrogate.regrfun = @regpoly0;
    surrogate.training = str2func([lower(surrogate.method) '_training']);
    surrogate.predictor = str2func([lower(surrogate.method) '_predictor']);
    surrogate.model = [];

metaproblem.par_objfun.surrogate = surrogate;

% function
metaproblem.objfun = @mask_objfun_mo_sa_outer; % depends on d and par_objfun. Needs to specify par_objfun.surrogate.model
return