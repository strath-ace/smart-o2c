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
problem.timeout_feasible = 120;

%% Generate first guess

x_guess = generate_guess(problem,1);

%% Solve problem

options = optimoptions(@fmincon,'Algorithm','sqp','Display','iter-detailed','MaxFunEvals',1e9,'MaxIter',1e9,'GradConstr','on','GradObj','on','TolFun',1e-6','TolCon',1e-6,'TolX',1e-12);
tstart = tic;
[x_sol,fval] = fmincon(@(x) multiphase_objectives(x,problem,1),x_guess,[],[],[],[],problem.norm_lb,problem.norm_ub,@(x) multiphase_constr_full(problem,x,1,tstart,-1),options);
toc

%% Plot solution

plot_solution_vs_time_multiphase(x_sol,problem,1,1)

[xt,ut,t] = eval_solution_over_time_multiphase(x_sol,problem,1);

xt = xt{1};
ut = ut{1};
t = t{1};
time = t;

%% Get individual variables to make plotting easier

alt = xt(:,1);
phi = xt(:,2);
theta = xt(:,3);
v = xt(:,4);
gamma = xt(:,5);
psi = xt(:,6);
r = alt+problem.structure{1}.constants.Re;

alpha = ut(:,1);
beta = ut(:,2);

ahat = alpha*180/pi;
bhat = beta*180/pi;
phihat = phi*180/pi;
thetahat = theta*180/pi;
gammahat = gamma*180/pi;
psihat = psi*180/pi;

rho = problem.structure{1}.constants.rho0*exp(-alt(:,1)/problem.structure{1}.constants.hr);

%% Other useful plots 

% thermal flux

qr =  17700*rho.^0.5.*(0.0001*v).^3.07;
qa = problem.structure{1}.constants.c0+problem.structure{1}.constants.c1*ahat+problem.structure{1}.constants.c2*ahat.^2+problem.structure{1}.constants.c3*ahat.^3;
figure ()
plot(time,qr.*qa)

% max gammadot

Cl = problem.structure{1}.constants.a0+problem.structure{1}.constants.a1*ahat;
L = 0.5*rho.*v.^2*problem.structure{1}.constants.S.*Cl;
grav = problem.structure{1}.constants.mu./(r.^2);

gammadot = L./(problem.structure{1}.constants.m.*v).*cos(beta) + cos(gamma).*(v./r-grav./v);           

figure()
plot(time,gammadot*180/pi)

% max psidot

psidot = L./(problem.structure{1}.constants.m.*v.*cos(gamma)).*sin(beta) + v./(r.*cos(theta)).*cos(gamma).*sin(psi).*sin(theta);

figure()
plot(time,psidot*180/pi)
