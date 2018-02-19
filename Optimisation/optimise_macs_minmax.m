% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [minmax, minmin] = optimise_macs_minmax(problem, algo_minmax, algo_outer, algo_inner)
% worst case scenario (min-max) and best case scenario (min-min): 
% - minimisation over the design variables d;
% - maximisation over the epistemic variables u;
% - constrained problem.
%
%% INPUT
%
% * problem: structure containing the informations of the problem:
%          * problem.output: choose to run minmax (0),  minmin (1) or both
%            minmax and minmin (2);
%          * problem.dim_u: dimention of the uncertain space;
%          * problem.lb_u{i}: lower boundaries of the uncertain variables
%            for the i-th objective function;
%          * problem.ub_u{i}: upper boundaries of the uncertain variables
%            for the i-th objective function;
%          * problem.dim_d: dimention of the design space;
%          * problem.lb_d: lower boundaries of the uncertain variables
%          * problem.ub_d: upper boundaries of the uncertain variables
%          * problem.fix: structure of the fixed parameters (if any);
%          * problem.n_obj: number of objective function;
%          * problem.objfun = {@...}: objective function;
%          * problem.constraints = {@...}: constraint function;
% * algo_minmax: structure containing MP-AIDEA specific input information
%                (if empty default values are used) for the meta-algorithm.
%              * algo_minmax.par_minmax.maxnfeval: max number of function
%                evaluation for the minmax problem;
% * algo_outer: structure containing MP-AIDEA specific input information
%                (if empty default values are used) for the minimisation loop.
%              * algo_outer.par_mpaidea.nFeValMax: max number of function 
%                evaluation for the outer loop (minimisation over d);
%              * algo_outer.par_mpaidea.n_populations;
%              * algo_outer.par_mpaidea.n_agents:
%              * algo_outer.par_mpaidea.max_LR: maximum number of local 
%                restart before global restart (for cases when only one 
%                population is considered and no adaptation of delta_local and
%                local/global restart is performed; if par_mpaidea.max_LR = []
%                 then adaptation is performed);
% * algo_inner: structure containing MP-AIDEA specific input information
%                (if empty default values are used) for the minimisation loop.
%             * algo_inner.par_mpaidea.nFeValMax;  
%             * algo_inner.par_mpaidea.n_populations;
%             * algo_inner.par_mpaidea.n_agents;
%             * algo_inner.par_mpaidea.max_LR;
%
%
%% OUTPUT
%
% * minmax: structure with the worst case solution:
%         * minmax.d;
%         * minmax.u;
%         * minmax.f;
% * minmin: structure with the best case solution:
%         * minmin.d;
%         * minmin.u;
%         * minmin.f;
%
%% Authors: Massimiliano Vasile, Carlos Ortegga Absil, Gianluca Filippi, Marilena Di Carlo
% email: 
%
%% References:


  
minmax = [];
minmin = [];

% run minmax and/or minmin
if problem.output == 0 || problem.output == 2
    [minmax] = evaluate_minmax(problem, algo_minmax, algo_outer, algo_inner);
end

if problem.output == 1 || problem.output == 2
    [minmin] = evaluate_minmin(problem, algo_minmax, algo_outer, algo_inner);
end




return

