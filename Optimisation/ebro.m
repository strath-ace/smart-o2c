% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [minmin, minmax, LIST, LIST_EXACT] = ebro(problem, algo_minmax, algo_outer, algo_inner, algo_decomposition)
% =========================================================================
% EBRO (Evidence-Based Robust Optimisation):
% - UQ using epistemic uncertainty;
% - Evidence theory with Belief and Plausibility curves;
% - Decomposition approach to treat with complex systems: uncertain 
%   variables that couple two or more sub-systems and uncertain variables 
%   that influence only one subsystems are studied separatedly.
% - The Decomposition approach assure conservative curves at  polynomial
%   cost under some assumptions (see the reference papers).
%
% Reference:
% M. Vasile, G. Filippi, C. O. Absil, A. Riccardi, "Fast Belief Estimation 
% in Evidence Network Models", EUROGEN 2017, September 13-15,Madrid, Spain.
% G. Filippi, M. Vasile, M. Marchi, P. Vercesi, "Evidence-Based Robust 
% Optimisation of Space Systems with Evidence Network Models", IEEE World
% Congress on Computational Intelligence, 8-13 July 2018, Rio de Janeiro,
% Brazil.
% =========================================================================
% =========================================================================
%
% INPUT
%
% * problem: structure containing the informations of the problem:
%          * problem.output: choose to reconstruct Belief (0),  Plausibility 
%            (1) or both Belief and Plausibility (2);
%          * problem.input: run minmax and/or minmin (0), load the design
%            vector d (1) or load the design vector d and the uncertain
%            vector u (2);
%          * problem.exact_curves: reconstruct Belief and/or Plausibility
%            curves with the decomposition approach but do not evaluate the
%            exact curve(s) (0), do bot decomposition and exact curve (1)
%            or do only exact curve(s) (2);
%          * problem.num_functions: number of sub-functions in which the
%            system can be decomposed;
%          * problem.num_samples: number of samples in the partial curves;
%            the bigger the more precise the decomposition curve; 
%          * problem.dim_u: dimention of the uncertain space [u1, u2, u12];
%          * problem.lb_u{i}: lower boundaries of the uncertain variables
%            for the i-th objective function;
%          * problem.ub_u{i}: upper boundaries of the uncertain variables
%            for the i-th objective function;
%          * problem.bpa: basic probability assignment of all the intervals
%            of each epistemic variable;
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
% * LIST: structur with the information of the Focal Elements considered in
%         the decomposition approach;
% * LIST EXACT: structur with the information of all the Focal Elements



[minmin, minmax, LIST, LIST_EXACT] = ebro_decomposition(problem, algo_minmax, algo_outer, algo_inner, algo_decomposition);


return
