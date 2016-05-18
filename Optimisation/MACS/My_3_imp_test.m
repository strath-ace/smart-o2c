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
    i
    tic
    [memory, nfeval,ener]=macs7v16(func,[],vlb,vub,opt,[],[]);
    toc
    tot = [tot;memory];
    %[M1(i),M2(i),M3(i),M4(i)]=ZDTmetrics(memory(:,1:end-4),memory(:,end-3:end-2),dom(:,1:end-4),dom(:,1:end-3:end-2),0,[1 1]);
    %[i M1(i) sum(M1<5)*100/i M4(i) sum(M4<5)*100/i]

%     i
%     figure(1000)
     %plot(memory(:,31),memory(:,32),'bo')
%     figure(2000)
%     plot(memory(:,3),memory(:,4),'bo')
     %axis equal
%     figure(3000)
%     plot(memory(:,3),memory(:,4),'bo')
%     plot(memory(:,4),memory(:,5),'bo')
%     axis equal
%     drawnow
%plot(memory(:,2),memory(:,3),'bo');
    
 end

