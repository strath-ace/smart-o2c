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

problem.timeout_feasible = 60;

%% MACS PARAMETERS

opt.maxnfeval=20e3;                                                         % maximum number of f evals
opt.popsize=10;                                                             % popsize (for each archive)
opt.rhoini=1;                                                               % initial span of each local hypercube (1=full domain)
opt.F=1;                                                                    % F, the parameter for Differential Evolution
opt.CR=1;                                                                   % CR, crossover probability
opt.p_social=1;                                                             % popratio
opt.max_arch=10;                                                            % archive size
opt.coord_ratio=1;
opt.contr_ratio=0.5;                                                        % contraction ratio
opt.draw_flag=0;                                                            % draw flag
opt.cp=0;                                                                   % constraints yes/no
opt.MBHflag=1;                                                              % number of MBH steps
opt.mbh_options = optimoptions(@fmincon,'Algorithm','sqp','Display','iter-detailed','MaxFunEvals',10*length(problem.norm_lb),'MaxIter',10*length(problem.norm_lb),'TolCon',1e-6,'TolFun',1e-6,'TolX',1e-9,'GradConstr','on','GradObj','on','DerivativeCheck','off');
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
opt.oc.problem = problem;

opt.timeout_single_level_opt = 3600;

tol_conv = 1e-6;
tol_fun = 1e3;
maxits = length(problem.norm_lb);
fminconoptions = optimoptions(@fmincon,'Algorithm','sqp','Display','off','MaxFunEvals',maxits,'MaxIter',maxits,'TolCon',tol_conv,'TolFun',tol_fun,'TolX',1e-12,'GradConstr','on','GradObj','on');

%% OPTIMISATION LOOP

fun = @(x_in) DFET_bilevel(x_in,1e7,[1e5;5e2],2,problem,fminconoptions);

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
    %[xt,ut,t] = eval_solution_over_time_multiphase(xx,problem,1);

%     for j = 1:problem.num_phases
%         
%         xxt = xt{j};
%         uut = ut{j};
%         tt = t{j};
%         time = tt;
%         
%         %% Get individual variables to make plotting easier
%         
%         x_sol2 = xx.*problem.scales.scale_opt;
%         static = x_sol2(logical((problem.other_vars).*(problem.phase_mask==j)));
%         S = static(1);
%         alt = xxt(:,1);
%         phi = xxt(:,2);
%         theta = xxt(:,3);
%         v = xxt(:,4);
%         gamma = xxt(:,5);
%         psi = xxt(:,6);
%         
%         w_e= problem.structure{1}.constants.omega_e;
%         
%         if j < 2
%             
%             m = xxt(:,7);
%             
%         else
%             
%             m = static(2);
%             
%         end
%         
%         r = alt+problem.structure{1}.constants.Re;
%         
%         alpha = uut(:,1);
%         beta = uut(:,2);
%         
%         ahat = alpha*180/pi;
%         bhat = beta*180/pi;
%         phihat = phi*180/pi;
%         thetahat = theta*180/pi;
%         gammahat = gamma*180/pi;
%         psihat = psi*180/pi;
%         
%         p = zeros(size(alt));
%         rho = zeros(size(alt));
%         c = zeros(size(alt));
%         Mach = zeros(size(alt));
%         Cl = zeros(size(alt));
%         Cd = zeros(size(alt));
%         E = zeros(size(alt));
%         
%         for k = 1:length(alt)
%             
%             [p(k),rho(k),c(k)] = atmo_ISA_smooth(alt(k));
%             Mach(k) = v(k)/c(k);
%             Cl(k) = problem.structure{j}.constants.cl_fun(ahat(k),Mach(k));
%             Cd(k) = problem.structure{j}.constants.cd_fun(ahat(k),Mach(k));
%             E(k) = Cl(k)./Cd(k);
%             
%         end
%         
%         if j < 2
%             
%             T = ( problem.structure{j}.constants.maxThrust-p* problem.structure{j}.constants.Ae*0.71082468).*uut(:,3);
%             
%         else
%             
%             T = 0;
%             
%         end
%         
%         
%         %% Other useful plots
%         
%         % % thermal flux
%         
%         v_imp = v*3.28084;      % convert velocity to imperial units
%         rho_imp = rho/515.379; % convert density to imperial units
%         
%         qr = 17700*rho_imp.^0.5.*(0.0001.*v_imp).^3.07;
%         qa = problem.structure{1}.constants.c0+problem.structure{1}.constants.c1*ahat+problem.structure{1}.constants.c2*ahat.^2+problem.structure{1}.constants.c3*ahat.^3;
%         q = qa.*qr.*11356.538527;     % convert heat flux density back into metric units
%         
%         % compute forces
%         
%         L = 0.5.*rho.*v.^2*S.*Cl;
%         D = 0.5.*rho.*v.^2*S.*Cd;
%         g = problem.structure{1}.constants.mu./(r.^2);
%         
%         % accelerations
%         
%         rdot = v.*cos(gamma);
%         vdot = (T.*cos(alpha)-D)./m - g.*sin(gamma)+w_e.^2.*r.*cos(theta).*(sin(gamma).*cos(theta)-cos(gamma).*sin(psi).*sin(theta));
%         gammadot = (T.*sin(alpha)+L)./(m.*v).*cos(beta) + cos(gamma).*(v./r-g./v)+2.*w_e.*cos(psi).*cos(theta)+w_e^2.*r./v.*cos(theta).*(sin(psi).*sin(gamma).*sin(theta)+cos(gamma).*cos(theta));
%         psidot = (T.*sin(alpha)+L)./(m.*v.*cos(gamma)).*sin(beta) + v./r.*cos(gamma).*cos(psi).*tan(theta)+2*w_e.*(sin(psi).*cos(theta).*tan(gamma)-sin(theta))-w_e^2.*r./(v.*cos(gamma)).*cos(theta).*sin(gamma).*cos(psi);
%         totacc2= vdot.^2+v.^2.*(gammadot.^2+psidot.^2);
%         
%         % absolute positions
%         xi = r.*cos(phi).*cos(theta);
%         yi = r.*sin(phi).*cos(theta);
%         zi = r.*sin(theta);
%         
%         figure(1)
%         plot(time,alt)
%         hold on
%         xlabel('t [s]')
%         ylabel('h [m/s]')        
%         
%         % vdot
%         
%         figure(2)
%         plot(time,vdot./9.81)
%         hold on
%         xlabel('t [s]')
%         ylabel('vdot [g]')
%         %plot([0 x_sol(problem.structure{1}.tf_vars)],[problem.structure{1}.constants.zeta*180/pi problem.structure{1}.constants.zeta*180/pi],'r')
%         %plot([0 x_sol(problem.structure{1}.tf_vars)],[-problem.structure{1}.constants.zeta*180/pi -problem.structure{1}.constants.zeta*180/pi],'r')
%         
%         % gammadot
%         
%         figure(3)
%         plot(time,v.*gammadot/9.81)
%         hold on
%         xlabel('t [s]')
%         ylabel('v*gammadot [g]')
%         %plot([0 x_sol(problem.structure{1}.tf_vars)],[problem.structure{1}.constants.zeta*180/pi problem.structure{1}.constants.zeta*180/pi],'r')
%         %plot([0 x_sol(problem.structure{1}.tf_vars)],[-problem.structure{1}.constants.zeta*180/pi -problem.structure{1}.constants.zeta*180/pi],'r')
%         
%         % psidot
%         
%         figure(4)
%         plot(time,v.*psidot/9.81)
%         hold on
%         xlabel('t [s]')
%         ylabel('v*psidot [g]')
%         %plot([0 x_sol(problem.structure{1}.tf_vars)],[problem.structure{1}.constants.zeta*180/pi problem.structure{1}.constants.zeta*180/pi],'r')
%         %plot([0 x_sol(problem.structure{1}.tf_vars)],[-problem.structure{1}.constants.zeta*180/pi -problem.structure{1}.constants.zeta*180/pi],'r')
%         
%         figure(5)
%         plot(time,totacc2.^0.5/9.81)
%         hold on
%         xlabel('t [s]')
%         ylabel('total acceleration [g]')
%         
%         figure(6)
%         plot(time,Mach)
%         hold on
%         xlabel('t [s]')
%         ylabel('Mach')
%         
%         figure(7)
%         plot(time,Cl./Cd)
%         hold on
%         xlabel('t [s]')
%         ylabel('Cl./Cd')
%         
%         figure(8)
%         plot(time,L,time,T.*sin(alpha),time,L+T.*sin(alpha))
%         hold on
%         xlabel('t [s]')
%         ylabel('Forces')
%         legend('L','T*sin(\alpha)','L+T*sin(\alpha)')
%         
%         figure (9)
%         plot(time,q)
%         hold on
%         xlabel('time [s]')
%         ylabel('q [W/m^2]')
%         
%         figure(10)
%         plot(problem.structure{j}.constants.Re*phi,alt)
%         hold on
%         xlabel('curvilinear abscissa [m]')
%         ylabel('altitude [m]')
%         
%         %     figure(6)
%         %     plot3(xi,yi,zi)
%         %     hold on
%         %     xlabel('x [m]')
%         %     ylabel('y [m]')
%         %     zlabel('z [m]')
%         
%     end
%     xx = xx.*structure.scale_optimisation_vars+structure.offset_optimisation_vars;
%        [x,u,xb] = extract_solution(xx(~structure.static_vars),structure,x_f);
%        t_f = xx(structure.tf_vars);
%        plot_solution_vs_time(x,u,x_0,xb,t_0,t_f,structure.uniform_els,structure,i+1)
%        subplot(2,1,1)
%        axis([0 t_f min(state_bounds(:,1)) max(state_bounds(:,2)) ])
%        subplot(2,1,2)
%        axis([0 t_f min(control_bounds(:,1)) max(control_bounds(:,2))])
%     drawnow
    
end

%% animation of Pareto front

% nit = max(mem.history(:,1));
% 
% mino1 = min(mem.history(:,end-3));
% maxo1 = max(mem.history(:,end-3));
% mino2 = min(mem.history(:,end-2));
% maxo2 = max(mem.history(:,end-2));
% deltao1 = maxo1-mino1;
% deltao2 = maxo2-mino2;
% 
% figure()
% 
% for i = 1:nit
% 
%     plot(mem.history(mem.history(:,1)==i,end-3),mem.history(mem.history(:,1)==i,end-2),'b*')
%     axis([mino1-deltao1*0.2 maxo1+deltao1*0.2 mino2-deltao2*0.2 maxo2+deltao2*0.2])
%     pause(0.5)
%     drawnow
% 
% end
% 
