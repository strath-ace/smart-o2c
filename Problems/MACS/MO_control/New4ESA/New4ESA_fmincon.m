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
problem.timeout_feasible = inf;

%% Generate first guess

x_guess = generate_guess(problem,1);

%x_guess = problem.norm_lb+rand(size(problem.norm_lb)).*(problem.norm_ub-problem.norm_lb);
%x_guess = (problem.norm_ub+problem.norm_lb)/2;
%x_guess(817) = (problem.norm_ub(317)-3*problem.norm_lb(317)/4);

% %% Solve problem
% 
% options = optimoptions(@fmincon,'Algorithm','interior-point','Display','iter-detailed','MaxFunEvals',1e4,'MaxIter',1e4,'GradConstr','on','GradObj','on','TolFun',1e3,'Tolcon',1e-6,'TolX',1e-12,'DerivativeCheck','off');%,'PlotFcns',@(x,optimValues,state) myoptimplot(x,optimValues,state,structure,t_0,x_0,x_f));
% tic
% x_sol = fmincon(@(x) feas_only(x),x_guess,[],[],[],[],problem.norm_lb,problem.norm_ub,@(x) multiphase_constr_full(problem,x,1),options);
% toc

options = optimoptions(@fmincon,'Algorithm','sqp','Display','iter-detailed','MaxFunEvals',1e9,'MaxIter',1e9,'GradConstr','on','GradObj','on','TolFun',1e-6','TolCon',1e-6,'TolX',1e-6,'DerivativeCheck','off');%,'PlotFcns',@(x,optimValues,state) myoptimplot(x,optimValues,state,structure,t_0,x_0,x_f));
tstart = tic;
[x_sol,fval,~,~,lambda,grad,hessian] = fmincon(@(x) multiphase_objectives(x,problem,1),x_guess,[],[],[],[],problem.norm_lb,problem.norm_ub,@(x) multiphase_constr_full(problem,x,1,tstart,-1),options);
%toc

%% Plot solution

plot_solution_vs_time_multiphase(x_sol,problem,1,1)

%qq = problem.structure{1}.uniform_in_nodes_state*t_fbest;
%time = sort(qq(:));
[xt,ut,t] = eval_solution_over_time_multiphase(x_sol,problem,1);

for j = 1:problem.num_phases
    
    xxt = xt{j};
    uut = ut{j};
    tt = t{j};
    time = tt;
    
    %% Get individual variables to make plotting easier
    
    x_sol2 = x_sol.*problem.scales.scale_opt;
    static = x_sol2(logical((problem.other_vars).*(problem.phase_mask==j)));

    S = static(1);    
    alt = xxt(:,1);
    phi = xxt(:,2);
    theta = xxt(:,3);
    v = xxt(:,4);
    gamma = xxt(:,5);
    psi = xxt(:,6);
    
    w_e= problem.structure{1}.constants.omega_e;
    
    if j < 3
        
        m = xxt(:,7);
        
    else
        
        m = static(2);
        
    end
    
    r = alt+problem.structure{1}.constants.Re;
    
    alpha = uut(:,1);
    beta = uut(:,2);
    
    if j<3
    
        delta = uut(:,3);
        
    else
        
        delta = 0;
        
    end
    
        
    ahat = alpha*180/pi;
    bhat = beta*180/pi;
    phihat = phi*180/pi;
    thetahat = theta*180/pi;
    gammahat = gamma*180/pi;
    psihat = psi*180/pi;
    
    p = zeros(size(alt));
    rho = zeros(size(alt));
    c = zeros(size(alt));
    Mach = zeros(size(alt));
    Cl = zeros(size(alt));
    Cd = zeros(size(alt));
    E = zeros(size(alt));
    
    for i = 1:length(alt)
        
        [p(i),rho(i),c(i)] = atmo_ISA_smooth(alt(i));
        Mach(i) = (v(i)-problem.structure{j}.constants.omega_e*r(i))/c(i);
        Cl(i) = problem.structure{j}.constants.cl_fun(ahat(i),Mach(i));
        Cd(i) = problem.structure{j}.constants.cd_fun(ahat(i),Mach(i));
        E(i) = Cl(i)./Cd(i);
        
    end
    
    T = ( problem.structure{j}.constants.maxThrust-p* problem.structure{j}.constants.Ae*0.71082468).*delta;
    
    %% Other useful plots
    
    % % thermal flux
    
    v_imp = (v-problem.structure{j}.constants.omega_e*r)*3.28084;      % convert velocity to imperial units
    rho_imp = rho/515.379; % convert density to imperial units

    qr = 17700*rho_imp.^0.5.*(0.0001.*v_imp).^3.07;
    qa = problem.structure{1}.constants.c0+problem.structure{1}.constants.c1*ahat+problem.structure{1}.constants.c2*ahat.^2+problem.structure{1}.constants.c3*ahat.^3;
    q = qa.*qr.*11356.538527;     % convert heat flux density back into metric units
    
    % compute forces
    
    L = 0.5.*rho.*(v-problem.structure{j}.constants.omega_e*r).^2*S.*Cl;
    D = 0.5.*rho.*(v-problem.structure{j}.constants.omega_e*r).^2*S.*Cd;
    g = problem.structure{1}.constants.mu./(r.^2);
        
    % accelerations

    rdot = v.*cos(gamma);
    
    vdot = (T.*cos(alpha)-D)./m - g.*sin(gamma);
    gammadot = (T.*sin(alpha)+L)./(m.*v).*cos(beta) + cos(gamma).*(v./r-g./v);
    psidot = (T.*sin(alpha)+L)./(m.*v.*cos(gamma)).*sin(beta) + v./(r.*cos(theta)).*cos(gamma).*sin(psi).*sin(theta);
     
    totacc2= vdot.^2+v.^2.*(gammadot.^2+psidot.^2);

    % absolute positions
    xi = r.*cos(theta).*cos(phi);
    yi = r.*cos(theta).*sin(phi);
    zi = r.*sin(theta);    
    
    % vdot
    
    figure(2)
    plot(time,vdot./9.81)
    hold on
    xlabel('t [s]')
    ylabel('vdot [g]')
    %plot([0 x_sol(problem.structure{1}.tf_vars)],[problem.structure{1}.constants.zeta*180/pi problem.structure{1}.constants.zeta*180/pi],'r')
    %plot([0 x_sol(problem.structure{1}.tf_vars)],[-problem.structure{1}.constants.zeta*180/pi -problem.structure{1}.constants.zeta*180/pi],'r')
    
    % gammadot
    
    figure(3)
    plot(time,v.*gammadot/9.81)
    hold on
    xlabel('t [s]')
    ylabel('v*gammadot [g]')    
    %plot([0 x_sol(problem.structure{1}.tf_vars)],[problem.structure{1}.constants.zeta*180/pi problem.structure{1}.constants.zeta*180/pi],'r')
    %plot([0 x_sol(problem.structure{1}.tf_vars)],[-problem.structure{1}.constants.zeta*180/pi -problem.structure{1}.constants.zeta*180/pi],'r')
    
    % psidot
    
    figure(4)
    plot(time,v.*psidot/9.81)
    hold on
    xlabel('t [s]')
    ylabel('v*psidot [g]')
    %plot([0 x_sol(problem.structure{1}.tf_vars)],[problem.structure{1}.constants.zeta*180/pi problem.structure{1}.constants.zeta*180/pi],'r')
    %plot([0 x_sol(problem.structure{1}.tf_vars)],[-problem.structure{1}.constants.zeta*180/pi -problem.structure{1}.constants.zeta*180/pi],'r')
    
    figure(5)
    plot(time,totacc2.^0.5/9.81)
    hold on
    xlabel('t [s]')
    ylabel('total acceleration [g]')
    
    figure(6)
    plot(time,Mach)
    hold on
    xlabel('t [s]')
    ylabel('Mach')
    
    figure(7)
    plot(time,Cl./Cd)
    hold on
    xlabel('t [s]')
    ylabel('Cl./Cd')
        
    figure(8)
    plot(time,gamma*180/pi)
    hold on
    xlabel('t [s]')
    ylabel('\gamma [deg]')    
    
    figure (9)
    plot(time,q)
    hold on
    xlabel('time [s]')
    ylabel('q [W/m^2]')
    
    figure(10)
    plot(problem.structure{j}.constants.Re*phi,alt)
    hold on
    xlabel('curvilinear abscissa [m]')
    ylabel('altitude [m]')
    
    figure(11)
    plot3(xi,yi,zi)
    hold on
    xlabel('x [m]')
    ylabel('y [m]')
    zlabel('z [m]')

end

% [xs,ys,zs] = sphere(100);
% xs = problem.structure{1}.constants.Re*xs;
% ys = problem.structure{1}.constants.Re*ys;
% zs = problem.structure{1}.constants.Re*zs;
% 
% figure(11)
% surf(xs,ys,zs)