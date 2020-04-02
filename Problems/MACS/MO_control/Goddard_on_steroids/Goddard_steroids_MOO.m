% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --a------------------
%

close all
clear
clc

%% Ensure calling script from it's folder

dirname = fileparts(mfilename('fullpath'));
cd(dirname);

%% Initialise problem

problem = initialise_problem(dirname);
problem.timeout_feasible = -1;

%% MACS PARAMETERS

opt.maxnfeval=20e3;                                                         % maximum number of f evals 
opt.popsize=20;                                                             % popsize (for each archive)
opt.rhoini=1;                                                               % initial span of each local hypercube (1=full domain)
opt.F=1;                                                                    % F, the parameter for Differential Evolution
opt.CR=1;                                                                   % CR, crossover probability
opt.p_social=1;                                                             % popratio
opt.max_arch=20;                                                            % archive size
opt.coord_ratio=1;                                              
opt.contr_ratio=0.5;                                                        % contraction ratio
opt.draw_flag=0;                                                            % draw flag
opt.cp=0;                                                                   % constraints yes/no 
opt.MBHflag=1;                                                              % number of MBH steps
opt.mbh_options = optimoptions(@fmincon,'Algorithm','sqp','Display','notify-detailed','MaxFunEvals',1000*length(problem.norm_lb),'MaxIter',1000*length(problem.norm_lb),'TolCon',1e-6,'TolFun',1e-6,'TolX',1e-9,'GradConstr','on','GradObj','on','DerivativeCheck','off');
opt.refine_freq = 10;
opt.smooth_scal_constr_fun = @Pascoletti_Serafini_constraints;
opt.cpat=0;                                                                 % pattern to DE
opt.explore_DE_strategy = 'rand';                                           % DE for exploring agents should pull towards the element with the best objective function value
opt.social_DE_strategy ='DE/current-to-rand/1';                             % DE for social agents
opt.explore_all = 1;                                                        % all agents should perform local search
opt.v = 0;
opt.int_arch_mult=1;
opt.dyn_pat_search = 1;
opt.upd_subproblems = 0;
opt.max_rho_contr = 5;
opt.pat_search_strategy = 'standard';
opt.vars_to_opt = problem.static_vars+problem.control_vars;         
opt.optimal_control = 1;
opt.bilevel = 1;
% optimal control settings (needed by DFET, included as structure of
% parameters of MACS just to make it simpler to pass and access them)
opt.oc.keepstatic = 0;
opt.oc.problem = problem;

tol_conv = 1e-7;
tol_fun = 1e-3;
maxits = 10*length(problem.norm_lb);
fminconoptions = optimoptions(@fmincon,'Algorithm','sqp','Display','off','MaxFunEvals',maxits,'MaxIter',maxits,'TolCon',tol_conv,'TolFun',tol_fun,'TolX',1e-12,'GradConstr','on','GradObj','off');

%% OPTIMISATION LOOP

fun = @(x_in) DFET_bilevel(x_in,1e7,[problem.structure{1}.tf_bounds(2);problem.structure{1}.xf_bounds(2,2)],1,problem,fminconoptions);

vlb = problem.lb;
vub = problem.ub;

for i=1:1
    
    [mem(i).memory,~,~,mem(i).history]=macs7v16OC(fun,[],vlb,vub,opt,[],[]);
    
end


%% plot 

[~,b] = sort(mem(1).memory(:,length(vlb)+2));

qq = mem(1).memory(b,:);    %sort wrt t_f

for i = 1:size(qq,1)
    
    xx = qq(i,1:length(vlb));
    xx = xx(:);
    xx = xx./problem.scales.scale_opt;          % macs stores dimensional solutions
    plot_solution_vs_time_multiphase(xx,problem,1,i+1)
    drawnow
    
end
