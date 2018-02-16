% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%
close all
clear
clc

%% Ensure calling script from it's folder

dirname = fileparts(mfilename('fullpath'));
cd(dirname);

%% Initialise problem

problem = initialise_problem(dirname);

%% Generate first guess

x_guess = generate_guess(problem,1);
problem.timeout_feasible = inf;

%% Solve problem

options = optimoptions(@fmincon,'Algorithm','sqp','Display','iter-detailed','MaxFunEvals',1e9,'MaxIter',1e9,'GradConstr','on','GradObj','on','TolFun',1e-6','TolCon',1e-6,'TolX',1e-12,'DerivativeCheck','off');%,'PlotFcns',@(x,optimValues,state) myoptimplot(x,optimValues,state,structure,t_0,x_0,x_f));
tstart = tic;
[x_sol,fval] = fmincon(@(x) multiphase_objectives(x,problem,1),x_guess,[],[],[],[],problem.norm_lb,problem.norm_ub,@(x) multiphase_constr_full(problem,x,1,tstart,-1),options);
toc

%% Plot solution

plot_solution_vs_time_multiphase(x_sol,problem,1)

[xt,ut,t] = eval_solution_over_time_multiphase(x_sol,problem,1);