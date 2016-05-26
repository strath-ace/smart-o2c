%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example of run of optimisation problem of CEC 2009 using MACS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

% Add path to optimiser folder
addpath(genpath('..\..\Optimisation'))

% Add path to % Add path to problem folder
addpath(genpath('CEC2009'))
tic

test_function = 'UF1';%'UF2';'UF3';'UF4';'UF5';'UF6';'UF7';'UF8';'UF9';'UF10'};

switch test_function
    case 'UF1'
        problem='cec09';
        np=30;
        no=2;
        name='UF1';
        vlb0=zeros(1,np);
        vub0=ones(1,np);
%         load CEC2009/UF1.dat
%         fp=UF1;
%         xp=fp;
%         weights=ones(1,2);
    case 'UF2'
        problem='cec09';
        np=30;
        no=2;
        name='UF2';
        vlb0=zeros(1,np);
        vub0=ones(1,np);
%         load CEC2009/UF2.dat
%         fp=UF2;
%         xp=fp;
%         weights=ones(1,2);
    case 'UF3'
        problem='cec09';
        np=30;
        no=2;
        name='UF3';
        vlb0=zeros(1,np);
        vub0=ones(1,np);
%         load CEC2009/UF3.dat
%         fp=UF3;
%         xp=fp;
%         weights=ones(1,2);
    case 'UF4'
        problem='cec09';
        np=30;
        no=2;
        name='UF4';
        vlb0=zeros(1,np);
        vub0=ones(1,np);
%         load CEC2009/UF4.dat
%         fp=UF4;
%         xp=fp;
%         weights=ones(1,2);
    case 'UF5'
        problem='cec09';
        np=30;
        no=2;
        name='UF5';
        vlb0=zeros(1,np);
        vub0=ones(1,np);
%         load CEC2009/UF5.dat
%         fp=UF5;
%         xp=fp;
%         weights=ones(1,2);
    case 'UF6'
        problem='cec09';
        np=30;
        no=2;
        name='UF6';
        vlb0=zeros(1,np);
        vub0=ones(1,np);
%         load CEC2009/UF6.dat
%         fp=UF6;
%         xp=fp;
%         weights=ones(1,2);
    case 'UF7'
        problem='cec09';
        np=30;
        no=2;
        name='UF7';
        vlb0=zeros(1,np);
        vub0=ones(1,np);
%         load CEC2009/UF7.dat
%         fp=UF7;
%         xp=fp;
%         weights=ones(1,2);
    case 'UF8'
        problem='cec09';
        np=30;
        no=3;
        name='UF8';
        vlb0=zeros(1,np);
        vub0=ones(1,np);
%         load CEC2009/UF8.dat
%         fp=UF8;
%         xp=fp;
%         weights=ones(1,3);
    case 'UF9'
        problem='cec09';
        np=30;
        no=3;
        name='UF9';
        vlb0=zeros(1,np);
        vub0=ones(1,np);
%         load CEC2009/UF9.dat
%         fp=UF9;
%         xp=fp;
%         weights=ones(1,3);
    case 'UF10'
        problem='cec09';
        np=30;
        no=3;
        name='UF10';
        vlb0=zeros(1,np);
        vub0=ones(1,np);
%         load CEC2009/UF10.dat
%         fp=UF10;
%         xp=fp;
%         weights=ones(1,3);
end

%% MACS2 SETTINGS

opt.maxnfeval=300000;     % maximum number of f evals
opt.popsize=150;         % popsize (forM4 each archive)
opt.rhoini=1;           % initial span of each local hypercube (1=full domain)
opt.F=0.9;              % F, the parameter for Differential Evolution
opt.CR=0.9;             % CR, crossover probability
opt.p_social=0.2;         % popratio
opt.max_arch = 100;

if no==3
    
    opt.max_arch=150;          % archive size
    
end

opt.coord_ratio=1;
opt.contr_ratio=0.5;         % contraction ratio
opt.draw_flag=0;          % draw flag
opt.cp=0;          % constraints yes/no
opt.MBHflag=0;          % number of MBH steps
opt.explore_DE_strategy = 'rand';
opt.social_DE_strategy = 'DE/current-to-rand/1';
opt.v = 0;
opt.dyn_pat_search = 1;
opt.upd_subproblems = 0;
opt.max_rho_contr = 5;
opt.pat_search_strategy = 'standard';
opt.optimal_control = 0;
opt.vars_to_opt = ones(length(vlb0),1);

[x,fval,exitflag,output] = optimise_macs(@test_func_cec09,vlb0,vub0,opt,name,np);

toc