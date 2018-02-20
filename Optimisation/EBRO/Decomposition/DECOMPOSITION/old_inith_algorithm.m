% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [algo_inner] = init_algorithm(problem)
[~, ~, algo_inner] = init_algo_so(problem);

algo_inner.par.nFeValMax = 400*problem.dim_u;            % number of function evaluations for IDEA

% %% ALGORITHM INNER LOOP: IDEA2
% algo_inner.optimise = @optimise_idea2_wrapper;                % algorithm in the form [x,f,exitflag,output] = algo(problem,par_algo)
%     %parameters{
%         options_idea2(1) = 300*problem.dim;                   % number of function evaluations for IDEA
%         options_idea2(4) = 5;                                 % number of agents
%         options_idea2(31) = 1;                                % number of local restarts
%         options_idea2(32) = 1;                                % communication heuristics (see com_red.m)
%         options_idea2(22) = 0.1;                              % delta_c: crowding factor for global restart
%         options_idea2(27) = 0.25;                             % tolconv: local convergence
%         options_idea2(15) = -1;                               % verbousity
%         options_idea2(28) = 0.2;                              % F
%         options_idea2(29) = 0.8;                              % CR
%         % options_idea2(30) = minmax_problem.n_obj;           % number of objectives (not used)
%         options_idea2(33) = 0.2;                              % initial probability of running the subproblem (not used)
% 
% algo_inner.par.options = options_idea2;
% algo_inner.par.initial_population = [];

end