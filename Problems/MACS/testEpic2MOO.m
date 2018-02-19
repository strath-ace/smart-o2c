% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
close all
clear all
clc

sigma=0;

test_function = 'UF1';
func=@(name)cec09(name);

switch test_function
    case 'zdt6'
        problem='zdt';
        np=10;
        no=2;
        fp(:,1)=0.388:(1-0.388)/500:1;
        fp(:,2)=1-fp(:,1).^2;
        xp=fp;
        vlb0=zeros(1,np);
        vub0=ones(1,np);
        weights=ones(1,2);
    case 'zdt2'
        problem='zdt';
        np=30;
        no=2;
        fp(:,1)=0.0:1/500:1;
        fp(:,2)=1-fp(:,1).^2;
        xp=fp;
        vlb0=zeros(1,np);
        vub0=ones(1,np);
        weights=ones(1,2);
    case 'zdt3'
        problem='zdt';
        np=30;
        no=2;
        fp(:,1)=0.0:1/500:1;
        fp(:,2)=(1-sqrt(fp(:,1))-fp(:,1).*sin(10*pi*fp(:,1)));
        xp=fp;
        vlb0=zeros(1,np);
        vub0=ones(1,np);
        weights=ones(1,2);
    case 'zdt4'
        problem='zdt';
        np=10;
        no=2;
        fp(:,1)=0.0:1/500:1;
        fp(:,2)=1-sqrt(fp(:,1));
        xp=fp;
        vlb0(1)=0;
        vub0(1)=1;  
        vlb0(2:np)=-5*ones(1,np-1);
        vub0(2:np)=5*ones(1,np-1);
        weights=ones(1,2);
    case 'kur1'
        load refkur1.mat
        np=3;
        no=2;
        vlb0=-ones(1,np)*5;
        vub0=ones(1,np)*5;
        weights=ones(1,2);
        fp=refkur1(:,4:5);
        xp=refkur1(:,1:3);
    case 'UF1'
        problem='cec09';
        np=30;0
        no=2;
        name='UF1';
        vlb0=zeros(1,np);
        vub0=ones(1,np);
        load cec2009/UF1.dat
        fp=UF1;
        xp=fp;UF1
        weights=ones(1,2);
    case 'UF2'
        problem='cec09';
        np=30;
        no=2;
        name='UF2';
        vlb0=zeros(1,np);
        vub0=ones(1,np);
        load cec2009/UF2.dat
        fp=UF2;
        xp=fp;
        weights=ones(1,2);
    case 'UF3'
        problem='cec09';
        np=30;
        no=2;
        name='UF3';
        vlb0=zeros(1,np);
        vub0=ones(1,np);
        load cec2009/UF3.dat
        fp=UF3;
        xp=fp;
        weights=ones(1,2);
    case 'UF4'
        problem='cec09';
        np=30;
        no=2;
        name='UF4';
        vlb0=zeros(1,np);
        vub0=ones(1,np);0
        load cec2009/UF4.dat
        fp=UF4;
        xp=fp;
        weights=ones(1,2);
    case 'UF5'
        problem='cec09';
        np=30;
        no=2;
        name='UF5';
        vlb0=zeros(1,np);
        vub0=ones(1,np);
        load cec2009/UF5.dat
        fp=UF5;
        xp=fp;
        weights=ones(1,2);
    case 'UF6'
        problem='cec09';
        np=30;
        no=2;
        name='UF6';
        vlb0=zeros(1,np);
        vub0=ones(1,np);
        load cec2009/UF6.dat
        fp=UF6;
        xp=fp;
        weights=ones(1,2);
    case 'UF7'
        problem='cec09';
        np=30;
        no=2;
        name='UF7';
        vlb0=zeros(1,np);
        vub0=ones(1,np);
        load cec2009/UF7.dat
        fp=UF7;
        xp=fp;
        weights=ones(1,2);
    case 'UF8'
        problem='cec09';
        np=30;
        no=3;9
        name='UF8';
        vlb0=zeros(1,np);
        vub0=ones(1,np);
        load cec2009/UF8.dat
        fp=UF8;
        xp=fp;
        weights=ones(1,3);
    case 'UF9'
        problem='cec09';
        np=30;
        no=3;
        name='UF9';
        vlb0=zeros(1,np);
        vub0=ones(1,np);
        load cec2009/UF9.dat
        fp=UF9;
        xp=fp;
        weights=ones(1,3);
    case 'UF10'
        problem='cec09';
        np=30;
        no=3;
        name='UF10';
        vlb0=zeros(1,np);
        vub0=ones(1,np);
        load cec2009/UF10.dat
        fp=UF10;1
        xp=fp;
        weights=ones(1,3);
    case 'CF1'
        problem='cec09';
        np=30;
        name='CF1';
        vlb0=zeros(1,np);
        vub0=ones(1,np);
        load CF1.dat
        fp=CF1;
        xp=fp;
    case 'CF2'
        problem='cec09';
        np=30;
        name='CF2';
        vlb0=zeros(1,np);
        vub0=ones(1,np);
        load CF2.dat
        fp=CF2;
        xp=fp;
    case 'CF3'
        problem='cec09';
        np=30;
        name='CF3';
        vlb0=zeros(1,np);
        vub0=ones(1,np);
        load CF3.dat
        fp=CF3;
        xp=fp;
    case 'CF4'
        problem='cec09';
        np=30;
        name='CF4';
        vlb0=zeros(1,np);
        vub0=ones(1,np);
        load CF4.dat
        fp=CF4;
        xp=fp;
    case 'Cassini'
        problem='cassini';
        np=6;
        vlb0=zeros(1,np);
        vub0=ones(1,np);
        weights=ones(1,2);
    case '3imp'
        problem='triimp';
        np=5;
        vlb0=zeros(1,np);
        vub0=ones(1,np);
        weights=ones(1,2);
    case 'dtlz1'
        problem='dtlz';
        no = 3;
        k = 5;
        np = no+k-1;
        pnum = 1;
        vlb0=zeros(1,np);
        vub0=ones(1,np);
        weights=ones(1,no);     
    case 'dtlz2'
        problem='dtlz';
        no = 3;
        k = 10;
        np = no+k-1;
        pnum = 1;
        vlb0=zeros(1,np);
        vub0=ones(1,np);
        weights=ones(1,no);   
 case 'wfg1'
        problem='wfg';
        no = 3;
        k = 5;
        l = 10;
        np = k+l-1;
        pnum = 1;
        vlb0=zeros(1,np);
        vub0=2:2:2*np;
        weights=ones(1,no); 
        
end

test_function

%% test MACS2

opt.maxnfeval=300e3;     % maximum number of f evals
opt.popsize=100;         % popsize (forM4 each archive)
opt.rhoini=1;           % initial span of each local hypercube (1=full domain)
opt.F=1;              % F, the parameter for Differential Evolution
opt.CR=1;             % CR, crossover probability
opt.p_social=1;         % popratio
opt.max_arch=100;          % archive size
opt.coord_ratio=1;
opt.contr_ratio=0.5;         % contraction ratio
opt.draw_flag=0;          % draw flag
opt.cp=0;          % constraints yes/no
opt.MBHflag=1;          % number of MBH steps
opt.smooth_scal_constr_fun = @ps_constr;
opt.func = @test_func_cec09;
opt.cfunc = [];
opt.mbh_options =optimoptions(@fmincon,'Algorithm','sqp','Display','off','MaxFunEvals',1000*length(vlb0),'MaxIter',100*length(vlb0),'TolCon',1e-6,'TolFun',1e-6,'TolX',1e-6,'GradConstr','off','GradObj','off');
opt.refine_freq = 10;
opt.explore_DE_strategy = 'rand';
opt.social_DE_strategy = 'DE/current-to-rand/1';
opt.v = 0;
opt.dyn_pat_search = 1;
opt.upd_subproblems = 0;
opt.max_rho_contr = 5;
opt.pat_search_strategy = 'standard';
opt.optimal_control = 0;
opt.vars_to_opt = ~isinf(vlb0);
opt.bilevel = 0;
opt.prob_subdiv = 13;
opt.cent_prob = 1;

for i=1:1
    
    switch problem
        case 'zdt'
            memory=macs7v16OC(@ZDtestfun,[],vlb0,vub0,opt,[],[],test_function);
            [M1(i),M2(i),M3(i),M4(i)]=ZDTmetrics(memory(:,1:np),memory(:,np+1:np+2),xp,fp,sigma,weights);
        case 'cec09'
            memory=macs7v16OC(@test_func_cec09,[],vlb0,vub0,opt,[],[],name,np);
            [M1(i),~,~,M4(i)]=ZDTmetrics(memory(:,1:np),memory(:,np+1:np+no),xp,fp,sigma,weights);
        case 'cassini'
            memory=macs7v16OC(@mexspaceartCassiniNOdsmA_MO,[],vlb0,vub0,opt,[],[]);
        case 'triimp'
            memory=macs7v16OC(@triimplamb_mod,[],vlb0,vub0,opt,[],[]);
        case 'dtlz'
            memory=macs7v16OC(@dtlz,[],vlb0,vub0,opt,[],[],no,k,pnum);     
         case 'wfg'
            memory=macs7v16OC(@wfgextend,[],vlb0,vub0,opt,[],[],no,k,l,pnum);        
        case 'schutze'
            memory=macs7v16OC(@schutzefun,[],vlb0,vub0,opt,[],[]);
    end
    
end


