% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------

close all
clear all
clc

tic

sigma=0;
trials = 200;

tf = {'UF1';'UF2';'UF3';'UF4';'UF5';'UF6';'UF7';'UF8';'UF9';'UF10';'zdt2';'zdt4';'zdt6';'Cassini'};
tau_spr = [4.5e-3;5e-3;2e-2;3e-2;5e-2;3e-2;4.5e-3;1;1;1;1;1;1;1];
tau_con = [5e-3;5e-3;2e-2;3.5e-2;3e-2;3e-2;5e-3;1;1;1;1;1;1;1];
v = 0;

% Create a name for a subfolder
MainFolder = sprintf('CEC-TEST-Corrected-Alpha2Trim-Normalized_Velocity');

% Create the folder if it doesn't exist already.
if ~exist(MainFolder, 'dir')
    mkdir(MainFolder);
end


for k = 1:size(v,1);
    
    for j = 1:size(tf,1);
        
        %% CREATE A SUBFOLDER FOR THE FUNCTION
        
        test_function = tf{j};
        
        fprintf ('Test function : %s\n', test_function);
        
        % Function Subfolder
        FunSubFolder = test_function;
        
        % Create the folder if it doesn't exist already.
        FunPath = sprintf('%s/%s',MainFolder,FunSubFolder);
        if ~exist(FunPath, 'dir')
            mkdir(FunPath);
        end
        
        clear M1 M4 vlb0 vub0 arch_all xp fp
        %test_function='zdt6';
        
        switch test_function
            case 'schwefel'
                problem='soo';
                np=10;
                no=1;
                fp(1)=0;
                xp=zeros(1,np);
                vlb0=zeros(1,np);
                vub0=ones(1,np);
                weights=ones(1,2);
            case 'rastrigin'
                problem='soo';
                np=10;
                no=1;
                fp(1)=0;
                xp=zeros(1,np);
                vlb0=zeros(1,np);
                vub0=ones(1,np);
                weights=ones(1,2);
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
                vlb0=-ones(1,np)*5;
                vub0=ones(1,np)*5;
                weights=ones(1,2);
                fp=refkur1(:,4:5);
                xp=refkur1(:,1:3);
            case 'UF1'
                problem='cec09';
                np=30;
                no=2;
                name='UF1';
                vlb0=zeros(1,np);
                vub0=ones(1,np);
                load cec2009/UF1.dat
                fp=UF1;
                xp=fp;
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
                vub0=ones(1,np);
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
                no=3;
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
                fp=UF10;
                xp=fp;
                weights=ones(1,3);
            case 'CF1'
                problem='cec09';
                np=30;
                name='CF1';
                vlb0=zeros(1,np);
                vub0=ones(1,np);
                load cec2009/CF1.dat
                fp=CF1;
                xp=fp;
            case 'CF2'
                problem='cec09';
                np=30;
                name='CF2';
                vlb0=zeros(1,np);
                vub0=ones(1,np);
                load cec2009/CF2.dat
                fp=CF2;
                xp=fp;
            case 'CF3'
                problem='cec09';
                np=30;
                name='CF3';
                vlb0=zeros(1,np);
                vub0=ones(1,np);
                load cec2009/CF3.dat
                fp=CF3;
                xp=fp;
            case 'CF4'
                problem='cec09';
                np=30;
                name='CF4';
                vlb0=zeros(1,np);
                vub0=ones(1,np);
                load cec2009/CF4.dat
                fp=CF4;
                xp=fp;
            case 'Cassini'
                problem='cassini';
                np=6;
                no=2;
                vlb0=zeros(1,np);
                vub0=ones(1,np);
                weights=ones(1,2);
            case '3imp'
                problem='triimp';
                np=5;
                no=2;
                vlb0=zeros(1,np);
                vub0=ones(1,np);
                weights=ones(1,2);
            case 'multi_const_test'
                problem='multi';
                np=2;
                vlb0=-ones(1,np);
                vub0=ones(1,np);
            case 'moo_robust'
                problem='robust';
                np=2;
                vlb0=zeros(1,np);
                vub0=ones(1,np)*10;
        end
        
        %% MACS2 SETTINGS
        
        if strcmp(problem,'cec09') || strcmp(problem,'cassini')
            
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
            opt.upd_subproblems = 1;
            opt.max_rho_contr = 5;
            opt.pat_search_strategy = 'standard';
            
        end
        
        if strcmp(problem,'zdt')
            
            opt.maxnfeval=15000;     % maximum number of f evals
            opt.popsize=10;         % popsize (forM4 each archive)
            opt.rhoini=1;           % initial span of each local hypercube (1=full domain)
            opt.F=0.9;              % F, the parameter for Differential Evolution
            opt.CR=0.9;             % CR, crossover probability
            opt.p_social=1;         % popratio
            opt.max_arch=200;          % archive size
            opt.coord_ratio=1;
            opt.contr_ratio=0.5;         % contraction ratio
            opt.draw_flag=0;          % draw flag
            opt.cp=0;          % constraints yes/no
            opt.MBHflag=0;          % number of MBH steps
            opt.explore_DE_strategy = 'rand';
            opt.social_DE_strategy = 'DE/current-to-rand/1';
            opt.v = 0;
            opt.dyn_pat_search = 1;
            opt.upd_subproblems = 1;
            opt.max_rho_contr = 5;
            opt.pat_search_strategy = 'standard';
            
        end
        
        %% CREATE A SUBFOLDER FOR THIS SETTING
        
        SetSubFolder = sprintf('%s',num2str(k));
        
        % Create the folder if it doesn't exist already.
        SetPath = sprintf('%s/Set-%s',FunPath,SetSubFolder);
        
        if ~exist(SetPath, 'dir')
            mkdir(SetPath);
        end
        
        %% RUN MANY INSTANCES
        
        % Save settings in a text file
        
        %// Extract field data
        fields = repmat(fieldnames(opt), numel(opt), 1);
        values = struct2cell(opt);
        
        %// Convert all numerical values to strings
        idx = cellfun(@isnumeric, values);
        values(idx) = cellfun(@num2str, values(idx), 'UniformOutput', 0);
        
        %// Combine field names and values in the same array
        C = {fields{:}; values{:}};
        fmt_str = repmat('%s = %s\n', 1, size(C, 2));
        
        filestring = sprintf('%s/settings.txt',SetPath);
        fileID = fopen(filestring,'w');
        fprintf(fileID,[fmt_str(1:end)], C{:});
        fclose(fileID);
        
        arch_all = [];
        
        M1 = zeros(trials,1);
        M2 = M1;
        M3 = M1;
        M4 = M1;
        n_agents_in_arch = M1;
        
        switch problem
            
            case 'cec09'
                
                parfor i=1:trials
                    
                    memories=macs7v16(@test_func_cec09,[],vlb0,vub0,opt,[],[],name,np);
                    
                    [M1(i),~,~,M4(i)]=ZDTmetrics(memories(:,1:np),memories(:,np+1:np+no),xp,fp,sigma,weights);
                    
                    memories = [i*ones(size(memories,1),1) memories  M1(i)*ones(size(memories,1),1) M4(i)*ones(size(memories,1),1)];
                    arch_all = [arch_all; memories];
                    n_agents_in_arch(i) = size(memories,1);
                    
                end
                
            case 'zdt'
                
                parfor i=1:trials
                    
                    memories=macs7v16(@ZDtestfun,[],vlb0,vub0,opt,[],[],test_function);
                    
                    [M1(i),~,~,M4(i)]=ZDTmetrics(memories(:,1:np),memories(:,np+1:np+no),xp,fp,sigma,weights);
                    
                    memories = [i*ones(size(memories,1),1) memories  M1(i)*ones(size(memories,1),1) M4(i)*ones(size(memories,1),1)];
                    arch_all = [arch_all; memories];
                    n_agents_in_arch(i) = size(memories,1);
                    
                end
                
            case 'cassini'
                
                parfor i=1:trials
                    
                    memories=macs7v16(@mexspaceartCassiniNOdsmA_MO,[],vlb0,vub0,opt,[],[]);
                    
                    %[M1(i),~,~,M4(i)]=ZDTmetrics(memories(:,1:np),memories(:,np+1:np+no),xp,fp,sigma,weights);
                    
                    memories = [i*ones(size(memories,1),1) memories  M1(i)*ones(size(memories,1),1) M4(i)*ones(size(memories,1),1)];
                    arch_all = [arch_all; memories];
                    n_agents_in_arch(i) = size(memories,1);
                    
                end
                
        end
        filestring = sprintf('%s/all-fronts.mat',SetPath);
        save(filestring,'arch_all');
        
        %% EVALUATE STATISTICS AND SCATTER PLOT
        
        % Unconditioned metrics
        mean_con = mean(M4);
        var_con = std(M4);
        mean_spr = mean(M1);
        var_spr = std(M1);
        S = sum(M1<tau_spr(j))*100/trials;
        C = sum(M4<tau_con(j))*100/trials;
        nfull = sum(n_agents_in_arch==opt.max_arch)*100/trials;
        
        % Unconditioned spread plot
        close all
        plot(M4,M1,'k.');
        title('Spread plot');
        xlabel('M4 metric: convergence error')
        ylabel('M1 metric: spread error')
        hold on
        line([tau_con(j) tau_con(j)],[0 tau_spr(j)],'Color','r');
        line([0 tau_con(j)],[tau_spr(j) tau_spr(j)],'Color','r');
        
        filestring = sprintf('%s/spread_plot.fig',SetPath);
        savefig(filestring)
        filestring = sprintf('%s/spread_plot.eps',SetPath);
        print(filestring,'-depsc')
        
        % Save metrics file
        filestring = sprintf('%s/metrics.txt',SetPath);
        fileID = fopen(filestring,'w');
        
        fprintf(fileID,'Unconditioned metrics:\n\n');
        fprintf(fileID,'Convergence: M4 mean = %s, M4 variance = %s\n',num2str(mean_con),num2str(var_con));
        fprintf(fileID,'Convergence probability below %s: %s\n', num2str(tau_con(j)), num2str(C));
        fprintf(fileID,'Spread: M1 mean = %s, M1 variance = %s\n',num2str(mean_spr),num2str(var_spr));
        fprintf(fileID,'Spread probability below %s : %s\n', num2str(tau_spr(j)), num2str(S));
        fprintf(fileID,'Percentage of full archives = %s\n\n',num2str(nfull));
        
        fclose(fileID);
        
    end
end

toc