% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
close all
clear
clc
% formulation from Betts

%% Definition of problem specific constants

constants.n = 0.06511*pi/180;    % orbital rotational rate [rad/s]
constants.hmax = 1e4;    % max angular momentum provided by cmg [ft*lbf*sec]
constants.umax = 200;
constants.J = [ 2.80701911616e7 4.822509936e5 -1.71675094448e7
                4.822509936e5 9.5144639344e7 6.02604448e4
                -1.71675094448e7 6.02604448e4 7.6594401336e7]; % ISS Inertia tensor [slug*ft^2]
constants.Jinv = inv(constants.J);  % avoid inverting it again, it's a constant!

constants.vcross = @(v) [0 -v(3) v(2);
                         v(3) 0 -v(1);
                         -v(2) v(1) 0];

% IC
r1i = 2.9963689649816e-3;
r2i = 1.5334477761054e-1;
r3i = 3.8359805613992e-3;

omega1i = -9.5380685844896e-6;
omega2i = -1.1363312657036e-3;
omega3i = 5.3472801108427e-6;

h1i = 5000;
h2i = 5000;
h3i = 5000;

% EC

h1f = 0;
h2f = 0;
h3f = 0;

%% Definition of bounds on states, controls, times (i.e. transcription variables)

state_bounds = [-1 1; -1 1; -1 1; -1e-2 1e-2; -1e-2 1e-2; -1e-2 1e-2; -1e4 1e4; -1e4 1e4; -1e4 1e4];     % q1, q2, q3, q4, w1, w2, w3, h1, h2, h3
control_bounds = [-200 200; -200 200; -200 200];               % ux, uy, uz
time_bounds = [0 1800];                        % this is the maximum range of time, t0 and tf must be within this range

%% Definition of initial conditions, final conditions, static variables, with their bounds if free

imposed_t_0 = 1;                                % t_0 is fixed
t_0 = 0;                                        % if t_0 is fixed, impose this value, otherwise it is ignored
t0_bounds = [];                                 % if empty, use time_bounds. Not used if t_0 is imposed

imposed_t_f = 1;                                % t_f is free
t_f = 1800;                                       % if t_f is fixed, impose this value, otherwise it is ignored
tf_bounds = [];                             % if empty, use time_bounds. Not used if t_f is imposed

imposed_initial_states = ones(9,1);           % mask vector with the variables on which an initial condition is imposed
x_0 = [r1i r2i r3i omega1i omega2i omega3i h1i h2i h3i]';                         
x0_bounds = [];                                 % if empty, use state_bounds

imposed_final_states = [0 0 0 0 0 0 1 1 1]';             % mask vector with the variables on which a final condition is imposed
x_f = [0 0 0 0 0 0 h1f h2f h3f]';                
xf_bounds = [];                

other_vars_bounds = [0 constants.hmax];       % bounds on gamma (max actuator effort)

%% DFET Discretisation settings

num_elems = 10;
state_order = 9;
control_order = 9;
DFET = 1;
state_distrib = 'Lobatto';
control_distrib = 'Lobatto';
num_eqs = size(state_bounds,1);
num_controls = size(control_bounds,1);

%% Checks (could be estended to check validity of bounds on states and i.c., free vs fixed i.c. etc)

test_order = state_order+(DFET==1);

if DFET==0

    total_num_equations = num_elems*(test_order+1)*num_eqs+sum(imposed_final_states);
    total_num_unknowns = num_elems*(state_order+1)*num_eqs+num_elems*(control_order+1)*num_controls;

else
   
    total_num_equations = ((test_order+1)+(test_order)*(num_elems-1))*num_eqs;
    total_num_unknowns = ((test_order)*num_elems)*num_eqs+num_elems*(control_order+1)*num_controls+sum(imposed_final_states==0);
    
end

if total_num_equations>total_num_unknowns
    
    error('Problem is over-constrained: either reduce number of final constraints, increase number of control variables, use higher order polynomials for control variables, or use more elements');
    
end

if total_num_equations==total_num_unknowns
    
    warning('Number of constraints equal to number of unknown control coefficients: optimal control is not possible, only constraints satisfaction');
    
end

if any(x_0(logical(imposed_initial_states))<state_bounds(logical(imposed_initial_states),1)) || any(x_0(logical(imposed_initial_states))>state_bounds(logical(imposed_initial_states),2))
    
    error('Imposed initial conditions are outside of bounds specified for states')
        
end

if any(x_f(logical(imposed_final_states))<state_bounds(logical(imposed_final_states),1)) || any(x_f(logical(imposed_final_states))>state_bounds(logical(imposed_final_states),2))
    
    error('Imposed final conditions are outside of bounds specified for states')
        
end

%% Generation of transcription (problem independent!)

structure = prepare_transcription(num_eqs,num_controls,num_elems,state_order,control_order,DFET,state_distrib,control_distrib);
structure = impose_boundary_conditions(structure,imposed_initial_states,imposed_final_states,imposed_t_0,imposed_t_f);

%% Compute bounds, gauge normalisation factors

structure = transcribe_bounds(x_0,x_f,t_0,t_f,state_bounds,control_bounds,time_bounds,t0_bounds,tf_bounds,x0_bounds,xf_bounds,other_vars_bounds,structure);

%% Definition of dynamics (can be moved after the initial transcription, so the scaling factors are computed internally and the user cannot mess them up, and they are ensured to be consistent through all the code)

% Wrapping commonly used constants into a handy and clean structure. 
% Could also include functions/models, like atmospheric models, etc
% Very useful and flexible.

f = @(x,u,t,static,scales) ZPM_alt_state_equations(x,u,t,static,scales,constants);
dfx = [];   % auto-compute :)
dfu = [];   % auto-compute :)

%% Definition of objective functions and Bolza's problem weights (as before)

g = @(x,u,t,static,scales) ZPM_alt_objective(x,u,t,static,scales,constants);
weights = [0 1; 1 0];
dgxf = [];   % auto-compute :)
dguf = [];   % auto-compute :)
dgxi = [];   % auto-compute :)
dgui = [];   % auto-compute :)

%% Definition of constraints

% Inequality path constraints

c = @(x,u,t,static,scales) ZPM_alt_ipath_constraints(x,u,t,static,scales,constants);
dcx = [];   % auto-compute :)
dcu = [];   % auto-compute :)

% Equality path constraints

e = [];
dex = [];   % auto-compute :)
deu = [];   % auto-compute :)

% Final time inequality constraints and integral inequality constraints

h = [];
dhxf = [];% auto-compute :)
dhuf = [];% auto-compute :)
dhxi = [];% auto-compute :)
dhui = [];% auto-compute :)
wh = [];

% Final time equality constraints and integral equality contraints

q = @(x,u,t,static,scales) ZPM_alt_eepath_constraints(x,u,t,static,scales,constants);;
dqxf = [];% auto-compute :)
dquf = [];% auto-compute :)
dqxi = [];% auto-compute :)
dqui = [];% auto-compute :)
wq = [ones(6,1) zeros(6,1)];

%% Include function handles, derivatives, objective functions, constraints and weights

structure = include_functions(structure,f,dfx,dfu,g,weights,dgxf,dguf,dgxi,dgui,c,dcx,dcu,e,dex,deu,h,wh,dhxf,dhuf,dhxi,dhui,q,wq,dqxf,dquf,dqxi,dqui);

%% Generate first guess

t0_guess = 0;
tf_guess = t_f;
static_guess = [other_vars_bounds(2)-1];

x0_guess = x_0;
control_guess = [0 0 0];                                       
u_nodes = repmat(reshape(ones((control_order+1),1)*control_guess,num_controls*(control_order+1),1),num_elems,1);

structure.init_type = 'DFET-all';

tic
x_guess = make_first_guess(x0_guess,t0_guess,tf_guess,static_guess,u_nodes,structure);
toc

% normalise all

if strcmp(structure.init_type,'CloneIC')

    x_guess = [static_guess;x_guess];

end
x_guess = x_guess./structure.scales.scale_opt;

% hard clipping of x_guess within limits

x_guess(x_guess<structure.norm_lbv) = structure.norm_lbv(x_guess<structure.norm_lbv);
x_guess(x_guess>structure.norm_ubv) = structure.norm_ubv(x_guess>structure.norm_ubv);
x_sol = x_guess;

%% Plots

% load data
load('mem.mat')

% sort in ascending order of |h|max
[~,b] = sort(mem(1).memory(:,end-2));
zz = mem(1).memory(b,:);
    
Colorset = varycolor(size(zz,1));

plot3d = 0;

if plot3d
    
    for i = 1:size(zz,1)
        
        % get data
        xx = zz(i,1:length(x_guess));
        t_fbest = xx(structure.tf_vars);
        x0_best = x_0;
        x0_best(~structure.imposed_initial_states) = xx(structure.x0_vars);
        xx = xx(~structure.static_vars);
        [x_best,u_best,x_b] = extract_solution(xx,structure,x_f);
        qq = structure.uniform_in_nodes_state*t_fbest;
        time = sort(qq(:));
        [xt,ut] = eval_solution_over_time(x_best,u_best,0,t_fbest,time,structure.uniform_els,structure);
        
        % get individual states over time
        alt = xt(:,1);
        phi = xt(:,2);
        theta = xt(:,3);
        v = xt(:,4);
        gamma = xt(:,5);
        psi = xt(:,6);
        r = alt+constants.Re;
        
        % get individual controls over time
        alpha = ut(:,1);
        beta = ut(:,2);
        
        ahat = alpha*180/pi;
        bhat = beta*180/pi;
        phihat = phi*180/pi;
        thetahat = theta*180/pi;
        gammahat = gamma*180/pi;
        psihat = psi*180/pi;
        
        % get density and thermal fluxes (path constraints)
        rho = constants.rho0*exp(-alt(:,1)/constants.hr);
        qr =  17700*rho.^0.5.*(0.0001*v).^3.07;
        qa = constants.c0+constants.c1*ahat+constants.c2*ahat.^2+constants.c3*ahat.^3;
        q = qr.*qa;
        
        Cl = constants.a0+constants.a1*ahat;
        Cd = constants.b0+constants.b1*ahat+constants.b2*ahat.^2;
        L = 0.5*rho.*v.^2*constants.S.*Cl;
        D = 0.5*rho.*v.^2*constants.S.*Cd;
        grav = constants.mu./(r.^2);
        
        gammadot = L./(constants.m.*v).*cos(beta) + cos(gamma).*(v./r-grav./v);
        psidot = L./(constants.m.*v.*cos(gamma)).*sin(beta) + v./(r.*cos(theta)).*cos(gamma).*sin(psi).*sin(theta);
        
        figure(1) % alt over t
        box on
        set(gca,'LooseInset',get(gca,'TightInset'));
        plot3(time,i*ones(size(time)),alt,'Color',Colorset(i,:));
        hold on
        xlabel('time [s]');
        ylabel('Solution id')
        zlabel('altitude [ft]');
        
        figure(2) % phi over t
        box on
        set(gca,'LooseInset',get(gca,'TightInset'));
        plot3(time,i*ones(size(time)), phi,'Color',Colorset(i,:));
        hold on
        xlabel('time [s]');
        ylabel('Solution id')
        zlabel('\phi [rad]');
        
        figure(3) % theta over t
        box on
        set(gca,'LooseInset',get(gca,'TightInset'));
        plot3(time,i*ones(size(time)), theta,'Color',Colorset(i,:));
        hold on
        xlabel('time [s]');
        ylabel('Solution id')
        zlabel('\theta [rad]');
        
        figure(4) % v over t
        box on
        set(gca,'LooseInset',get(gca,'TightInset'));
        plot3 (time,i*ones(size(time)), v,'Color',Colorset(i,:));
        hold on
        xlabel('time [s]');
        ylabel('Solution id')
        zlabel('v [ft/s]');
        
        figure(5) % gamma over t
        box on
        set(gca,'LooseInset',get(gca,'TightInset'));
        plot3 (time,i*ones(size(time)), gamma,'Color',Colorset(i,:));
        hold on
        xlabel('time [s]');
        ylabel('Solution id')
        zlabel('\gamma [rad]');
        
        figure(6) % psi over t
        box on
        set(gca,'LooseInset',get(gca,'TightInset'));
        plot3 (time,i*ones(size(time)), psi,'Color',Colorset(i,:));
        hold on
        xlabel('time [s]');
        ylabel('Solution id')
        zlabel('\psi [rad]');
        
        figure(7) % alpha over t
        box on
        set(gca,'LooseInset',get(gca,'TightInset'));
        plot3 (time,i*ones(size(time)), ahat,'Color',Colorset(i,:));
        hold on
        xlabel('time [s]');
        ylabel('Solution id')
        zlabel('\alpha [deg]');
        
        figure(8) % beta over t
        box on
        set(gca,'LooseInset',get(gca,'TightInset'));
        plot3 (time,i*ones(size(time)), bhat,'Color',Colorset(i,:));
        hold on
        xlabel('time [s]');
        ylabel('Solution id')
        zlabel('\beta [deg]');
        
        figure(9) % q over t
        box on
        set(gca,'LooseInset',get(gca,'TightInset'));
        plot3 (time,i*ones(size(time)), q,'Color',Colorset(i,:));
        hold on
        xlabel('time [s]');
        ylabel('Solution id')
        zlabel('q [btu]');
        
        figure(10) % gammadot over t
        box on
        set(gca,'LooseInset',get(gca,'TightInset'));
        plot3 (time,i*ones(size(time)), gammadot,'Color',Colorset(i,:));
        hold on
        xlabel('time [s]');
        ylabel('Solution id')
        zlabel('\gamma dot [rad/s]');
        
        figure(11) % psidot over t
        box on
        set(gca,'LooseInset',get(gca,'TightInset'));
        plot3 (time,i*ones(size(time)), psidot,'Color',Colorset(i,:));
        hold on
        xlabel('time [s]');
        ylabel('Solution id')
        zlabel('\psi dot [rad/s]');
        
        figure(12) % L/D over t
        box on
        set(gca,'LooseInset',get(gca,'TightInset'));
        plot3 (time,i*ones(size(time)), L./D,'Color',Colorset(i,:));
        hold on
        xlabel('time [s]');
        ylabel('Solution id')
        zlabel('L/D [rad/s]');

        figure(13) % 3D trajectories
        box on
        set(gca,'LooseInset',get(gca,'TightInset'));
        plot3 ((alt+constants.Re).*cos(psi).*cos(theta),(alt+constants.Re).*sin(psi).*cos(theta), (alt+constants.Re).*sin(theta),'Color',Colorset(i,:));
        hold on
        xlabel('x');
        ylabel('y')
        zlabel('z');
        
    end
    
    figure(13)
    str = {'Solution 1','Solution 2','Solution 3','Solution 4','Solution 5','Solution 6','Solution 7','Solution 8','Solution 9','Solution 10'};
    legend(str);
    
else
    
    for i = 1:size(zz,1)
        
        % get data
        xx = zz(i,1:length(x_guess));
        t_fbest = t_f;
        x0_best = x_0;
        x0_best(~structure.imposed_initial_states) = xx(structure.x0_vars);
        xx = xx(~structure.static_vars);
        [x_best,u_best,x_b] = extract_solution(xx,structure,x_f);
        qq = structure.uniform_in_nodes_state*t_fbest;
        time = sort(qq(:));
        [xt,ut] = eval_solution_over_time(x_best,u_best,0,t_fbest,time,structure.uniform_els,structure);
        
        % get individual states over time
        q1 = xt(:,1);
        q2 = xt(:,2);
        q3 = xt(:,3);
        omega1 = xt(:,4);
        omega2 = xt(:,5);
        omega3 = xt(:,6);
        h1 = xt(:,7);
        h2 = xt(:,8);
        h3 = xt(:,9);
                
        % get individual controls over time
        u1 = ut(:,1);
        u2 = ut(:,2);
        u3 = ut(:,3);
        
        hmax = (h1.^2+h2.^2+h3.^2).^0.5;
        uint = (u1.^2+u2.^2+u3.^2);
        
        figure(1) % q1 over t
        box on
        set(gca,'LooseInset',get(gca,'TightInset'));
        plot(time,q1,'Color',Colorset(i,:));
        hold on
        xlabel('time [s]');
        ylabel('q_1 [rad]');
        
        figure(2) % q1 over t
        box on
        set(gca,'LooseInset',get(gca,'TightInset'));
        plot(time,q2,'Color',Colorset(i,:));
        hold on
        xlabel('time [s]');
        ylabel('q_2 [rad]');
        
        figure(3) % q1 over t
        box on
        set(gca,'LooseInset',get(gca,'TightInset'));
        plot(time,q3,'Color',Colorset(i,:));
        hold on
        xlabel('time [s]');
        ylabel('q_3 [rad]');
        
        figure(4) % q1 over t
        box on
        set(gca,'LooseInset',get(gca,'TightInset'));
        plot(time,omega1,'Color',Colorset(i,:));
        hold on
        xlabel('time [s]');
        ylabel('\omega_1 [rad/s]');
        
        figure(5) % q1 over t
        box on
        set(gca,'LooseInset',get(gca,'TightInset'));
        plot(time,omega2,'Color',Colorset(i,:));
        hold on
        xlabel('time [s]');
        ylabel('\omega_2 [rad/s]');
        
        figure(6) % q1 over t
        box on
        set(gca,'LooseInset',get(gca,'TightInset'));
        plot(time,omega3,'Color',Colorset(i,:));
        hold on
        xlabel('time [s]');
        ylabel('\omega_3 [rad/s]');
       
        figure(7) % q1 over t
        box on
        set(gca,'LooseInset',get(gca,'TightInset'));
        plot(time,h1,'Color',Colorset(i,:));
        hold on
        xlabel('time [s]');
        ylabel('h_1 [lbf*ft*s]');

        
        figure(8) % q1 over t
        box on
        set(gca,'LooseInset',get(gca,'TightInset'));
        plot(time,h2,'Color',Colorset(i,:));
        hold on
        xlabel('time [s]');
        ylabel('h_2 [lbf*ft*s]');
        
        figure(9) % q1 over t
        box on
        set(gca,'LooseInset',get(gca,'TightInset'));
        plot(time,h3,'Color',Colorset(i,:));
        hold on
        xlabel('time [s]');
        ylabel('h_3 [lbf*ft*s]');
        
        figure(10) % q1 over t
        box on
        set(gca,'LooseInset',get(gca,'TightInset'));
        plot(time,u1,'Color',Colorset(i,:));
        hold on
        xlabel('time [s]');
        ylabel('u_1 [lbf*ft]');
        
        figure(11) % q1 over t
        box on
        set(gca,'LooseInset',get(gca,'TightInset'));
        plot(time,u2,'Color',Colorset(i,:));
        hold on
        xlabel('time [s]');
        ylabel('u_2 [lbf*ft*s]');
       
        figure(12) % q1 over t
        box on
        set(gca,'LooseInset',get(gca,'TightInset'));
        plot(time,u3,'Color',Colorset(i,:));
        hold on
        xlabel('time [s]');
        ylabel('u_3 [lbf*ft*s]');

        figure(13) % q1 over t
        box on
        set(gca,'LooseInset',get(gca,'TightInset'));
        plot(time,hmax,'Color',Colorset(i,:));
        hold on
        xlabel('time [s]');
        ylabel('|h| [lbf*ft*s]');

        figure(14) % q1 over t
        box on
        set(gca,'LooseInset',get(gca,'TightInset'));
        plot(time,uint,'Color',Colorset(i,:));
        hold on
        xlabel('time [s]');
        ylabel('u^Tu [lbf^2*ft^2*s^2]');

        
    end
    
    str = {'Solution 1','Solution 2','Solution 3','Solution 4','Solution 5','Solution 6','Solution 7','Solution 8','Solution 9','Solution 10'};
        
end