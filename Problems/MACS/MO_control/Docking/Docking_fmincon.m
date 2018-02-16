% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------

close all
clear
clc

%% Ensure calling script from it's folder

dirname = fileparts(mfilename('fullpath'));
cd(dirname);

%% Initialise problem

problem = initialise_problem(dirname);
problem.timeout_feasible = inf;

%% Generate first guess

x_guess = generate_guess(problem,1);

%% Solve problem

options = optimoptions(@fmincon,'Algorithm','sqp','Display','iter-detailed','MaxFunEvals',1e9,'MaxIter',1e9,'GradConstr','on','GradObj','on','TolFun',1e-6','TolCon',1e-6,'TolX',1e-12,'DerivativeCheck','off');%,'PlotFcns',@(x,optimValues,state) myoptimplot(x,optimValues,state,structure,t_0,x_0,x_f));
tstart = tic;
[x_sol,fval] = fmincon(@(x) multiphase_objectives(x,problem,1),x_guess,[],[],[],[],problem.norm_lb,problem.norm_ub,@(x) multiphase_constr_full(problem,x,1,tstart,-1),options);
toc

%% Plot solution

plot_solution_vs_time_multiphase(x_sol,problem,1)

[xt,ut,t] = eval_solution_over_time_multiphase(x_sol,problem,1);

%qq = problem.structure{1}.uniform_in_nodes_state*t_fbest;
time =t{1};
%time = linspace(0,t_fbest,100');
%[xt,ut] = eval_solution_over_time(x_best,u_best,0,t_fbest,time,structure.uniform_els,structure);

%% 3d plot

figure(2)

ra = zeros(length(time),3);
rb = zeros(length(time),3);

for i = 1:length(time)

    q = xt{1}(i,7:10);
    p = xt{1}(i,14:17);
    S = Q_fun(q(1),q(2),q(3),q(4))';
    T = Q_fun(p(1),p(2),p(3),p(4))';
    ra(i,:) = (S*problem.structure{1}.constants.avec)';
    rb(i,:) = (T*problem.structure{1}.constants.bvec)';

end

plot3(xt{1}(:,1),xt{1}(:,2),xt{1}(:,3),'b');
hold on
plot3(0,0,0,'ro');
axis equal
plot3(0,0,0,'r');
quiver3(xt{1}(:,1),xt{1}(:,2),xt{1}(:,3),ra(:,1),ra(:,2),ra(:,3),0);
quiver3(0,0,0,rb(end,1),rb(end,2),rb(end,3),0);

%% 3d anim


for i = 1:length(time)
    
    figure(100)
    plot3(xt{1}(1:i,1),xt{1}(1:i,2),xt{1}(1:i,3),'b');
    hold on
    plot3(xt{1}(i,1),xt{1}(i,2),xt{1}(i,3),'bo');
    plot3(0,0,0,'ro');
    quiver3(xt{1}(i,1),xt{1}(i,2),xt{1}(i,3),ra(i,1),ra(i,2),ra(i,3),0);
    quiver3(0,0,0,rb(i,1),rb(i,2),rb(i,3),0);
    axis([-10 10 -10 10 -10 10])
    hold off
    drawnow
    %pause(0.1)

end