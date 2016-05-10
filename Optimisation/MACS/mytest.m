close all
clear all
clc

opt.maxnfeval=15000;                                                       % maximum number of f evals 
opt.popsize=30;                                                             % popsize (for each archive)
opt.rhoini=1;                                                               % initial span of each local hypercube (1=full domain)
opt.F=1;                                                                  % F, the parameter for Differential Evolution
opt.CR=1;                                                                 % CR, crossover probability
opt.p_social=1;                                                             % popratio
opt.max_arch=200;                                                            % archive size
opt.coord_ratio=1;                                                          % not used
opt.contr_ratio=0.5;                                                        % contraction ratio
opt.draw_flag=0;                                                            % draw flag
opt.cp=0;                                                                   % constraints yes/no 
opt.MBHflag=0;                                                              % number of MBH steps
opt.cpat=0;                                                                 % pattern to DE
opt.explore_DE_strategy = 'rand';                                           % DE for exploring agents should pull towards the element with the best objective function value
opt.social_DE_strategy ='DE/current-to-rand/1';                             % DE for social agents
opt.explore_all = 1;                                                        % all agents should perform local search
opt.v = 0;
opt.int_arch_mult=1;

%func=@my_test_problem2;
%func=@cec2009test;
func=@ZDtestfun;

%narc = 2;

%func = @ose2_fun;
%vlb = 0*ones(1,1+4*narc);
%vlb(end-narc+1:end) = -1;
%vub = 1*ones(1,1+4*narc);

vlb = 0*ones(1,10);
vlb(1) = 0;
%vlb(1:2) = 0;
%vlb(3:4) = 1e-4;
%vlb(7:8) = 1e-4;
vub = 1*ones(1,10);
vub(1)= 1;

%func = @name()


% for i=1:50
    tic
    [memory, nfeval,ener]=macs7v16(func,[],vlb,vub,opt,[],[],'zdt6');
    toc
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
    
% end

