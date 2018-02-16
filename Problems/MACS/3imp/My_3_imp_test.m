% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%-----------Copyright (C) 2016 University of Strathclyde-------------

close all
clear all
clc

opt.maxnfeval=30000;                                                       % maximum number of f evals 
opt.popsize=20;                                                             % popsize (for each archive)
opt.rhoini=1;                                                               % initial span of each local hypercube (1=full domain)
opt.F=1;                                                                    % F, the parameter for Differential Evolution
opt.CR=1;                                                                   % CR, crossover probability
opt.p_social=1;                                                           % popratio
opt.max_arch=100;                                                            % archive size
opt.coord_ratio=1;                                                          % not used
opt.contr_ratio=0.5;                                                        % contraction ratio
opt.draw_flag=0;                                                            % draw flag
opt.cp=0;                                                                   % constraints yes/no 
opt.MBHflag=0;                                                              % number of MBH steps
opt.cpat=0;                                                                 % pattern to DE
opt.explore_DE_strategy = 'rand';                                           % DE for exploring agents should pull towards the element with the best objective function value
opt.social_DE_strategy ='DE/current-to-rand/1';                             % DE for social agents
opt.explore_all = 1;                                                        % all agents should perform local search
opt.v = 1;
opt.int_arch_mult=1;

func=@triimplamb_mod;

vlb = [0 0 0 0 0];
vub = [1 1 1 1 1];

tot = [];

%load('tri_imp_reference_pareto_front.mat');

for i=1:1
    
    tic
    [memory, nfeval,ener]=macs7v16OC(func,[],vlb,vub,opt,[],[]);
    toc
    tot = [tot;memory];
    
 end

