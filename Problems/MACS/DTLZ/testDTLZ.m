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

%% test MACS2

opt.rhoini=1;           % initial span of each local hypercube (1=full domain)
opt.F=0.9;              % F, the parameter for Differential Evolution
opt.CR=0.9;             % CR, crossover probability
opt.p_social=1;         % popratio
opt.coord_ratio=1;
opt.contr_ratio=0.5;         % contraction ratio
opt.draw_flag=0;          % draw flag
opt.cp=0;          % constraints yes/no
opt.MBHflag=0;          % number of MBH steps
opt.smooth_scal_constr_fun = @ps_constr;
opt.func = @dtlz;
opt.cfunc = [];
opt.refine_freq = 1;
opt.explore_DE_strategy = 'rand';
opt.social_DE_strategy = 'DE/current-to-rand/1';
opt.v = 0;
opt.dyn_pat_search = 1;
opt.upd_subproblems = 0;
opt.max_rho_contr = 5;
opt.pat_search_strategy = 'standard';
opt.optimal_control = 0;
opt.bilevel = 0;
opt.cent_prob = 1;

sd = [13,4,3,3];

for pnum=[5,7]

    for no = 3%[3,6,8]
        
        j = (no==3)+2*(no==6)+3*(no==8);
        
        k = 5*(pnum==1)+10*(pnum>1);        
        np = no+k-1;
        vlb0=zeros(1,np);
        vub0=ones(1,np); 
        opt.vars_to_opt = ~isinf(vlb0);
       
        opt.mbh_options =optimoptions(@fmincon,'Algorithm','sqp','Display','off','MaxFunEvals',100*length(vlb0),'MaxIter',100*length(vlb0),'TolCon',1e-6,'TolFun',1e-6,'TolX',1e-6,'GradConstr','off','GradObj','off');        
        opt.prob_subdiv = 2;%3;%13*(no==3)+4*(no==6)+3*(no>6);

        opt.popsize = size(simplex_lattice_design(no,opt.prob_subdiv),1)+1;
%         if no==3
%             opt.popsize = opt.popsize-1;
%         end
        opt.max_arch=105*(no==3)+132*(no==6)+156*(no==8)+275*(no==10);          % archive size
        opt.maxnfeval=opt.max_arch*1000/2*(mod(pnum,2)==1)+opt.max_arch*500/2*(mod(pnum,2)==0);     % maximum number of f evals
        
        for i=1:20
    
            memory{pnum,j,i}=macs7v16OC(@dtlz,[],vlb0,vub0,opt,[],[],no,k,pnum);
    
        end
        
    end

end


