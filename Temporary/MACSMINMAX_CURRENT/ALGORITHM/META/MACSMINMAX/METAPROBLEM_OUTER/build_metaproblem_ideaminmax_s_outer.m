function [metaproblem] = build_metaproblem_ideaminmax_s_outer(problem)

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

    if problem.n_obj == 1
        surrogate.indicator = str2func([lower(surrogate.method) '_EI2']);
    elseif problem.n_obj == 2
        surrogate.indicator = str2func([lower(surrogate.method) '_EIdom']);
    else
        error('We do not have kriging indicators for many-objective (yet)' );
    end


metaproblem.par_objfun.objectives = 1:problem.n_obj;
metaproblem.par_objfun.surrogate = surrogate;
metaproblem.par_objfun.ymin = [];
metaproblem.par_objfun.sign_inner = problem.sign_inner;
% function
metaproblem.objfun = @mask_objfun_ideaminmax_s_outer; % depends on d and par_objfun. Needs to specify par_objfun.surrogate.model, par_objfun.ymin
return