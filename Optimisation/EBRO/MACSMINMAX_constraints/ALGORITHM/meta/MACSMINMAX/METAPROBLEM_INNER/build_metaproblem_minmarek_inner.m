% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [metaproblem] = build_metaproblem_minmarek_inner(problem)

% chromosome
dim_u =  problem.dim_u;
metaproblem.dim = dim_u;
metaproblem.lb = zeros(1,dim_u);
metaproblem.ub = ones(1,dim_u);

% surrogate
% These are hardcoded (HC) at the moment but should really come from par_minmax for flexibility
    surrogate.method = 'kriging';
    surrogate.corrfun = @corrgauss;
    surrogate.regrfun = @regpoly0;
    surrogate.training = str2func([lower(surrogate.method) '_training']);
    surrogate.predictor = str2func([lower(surrogate.method) '_predictor']);
    surrogate.indicator = str2func([lower(surrogate.method) '_EI']);
    surrogate.model = [];
metaproblem.par_objfun.surrogate = surrogate;

metaproblem.par_objfun.sign = problem.sign_inner; %dunno if necessary
metaproblem.par_objfun.ymin = [];

metaproblem.objfun = @mask_objfun_minmarek_inner; %depends on u and par_objfun. Needs to specify par_objfun.surrogate.model, par_objfun.ymin

return