close all
clear
clc
% formulation from Betts

%% Definition of problem specific constants

constants.mu = 1;
constants.a = 1e-2;

% IC
ri = 1.1;
vri = 0;
thetai = 0;
vti = 1/sqrt(ri);

%% Definition of bounds on states, controls, times (i.e. transcription variables)

state_bounds = [0.1 400; -5 5; -pi 20*pi; -5 5];     % r, vr, theta, vt
control_bounds = [-pi pi];                       % u (rad)
time_bounds = [0 80];                        % this is the maximum range of time, t0 and tf must be within this range

%% Definition of initial conditions, final conditions, static variables, with their bounds if free

imposed_t_0 = 1;                                % t_0 is fixed
t_0 = 0;                                        % if t_0 is fixed, impose this value, otherwise it is ignored
t0_bounds = [];                                 % if empty, use time_bounds. Not used if t_0 is imposed

imposed_t_f = 0;                                % t_f is free
t_f = 40;                                       % if t_f is fixed, impose this value, otherwise it is ignored
tf_bounds = [20 80];                            % if empty, use time_bounds. Not used if t_f is imposed

imposed_initial_states = ones(4,1);           % mask vector with the variables on which an initial condition is imposed
x_0 = [ri vri thetai vti]';                         
x0_bounds = [];                                 % if empty, use state_bounds

imposed_final_states = zeros(4,1);             % mask vector with the variables on which a final condition is imposed
x_f = [0 0 0 0]';                
xf_bounds = [];                

other_vars_bounds = [];                         % no static vars

%% DFET Discretisation settings

num_elems = 30;
state_order = 1;
control_order = 1;
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

f = @(x,u,t,static,scales) orbit_state_equations(x,u,t,static,scales,constants);
dfx = [];   % auto-compute :)
dfu = [];   % auto-compute :)

%% Definition of objective functions and Bolza's problem weights (as before)

g = @(x,u,t,static,scales) orbit_objective(x,u,t,static,scales,constants);
weights = [1 0; 1 0];
dgxf = [];   % auto-compute :)
dguf = [];   % auto-compute :)
dgxi = [];   % auto-compute :)
dgui = [];   % auto-compute :)

%% Definition of constraints

% Inequality path constraints

c = [];
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

q = [];
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
static_guess = [];

x0_guess = x_0;
control_guess = [0];
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

%load data
%load('/home/la/Dropbox/smart-o2c/Problems/MACS/MO_control/Orbit/New_feas_New_Utopia/selection_e-7tol.mat')
%mem.memory = selection;

load('mem_50k_e-9tol.mat')

ar = 1;
% sort in ascending order of |h|max
[~,b] = sort(mem(ar).memory(:,end-2));
zz = mem(ar).memory(b,:);

Colorset = varycolor(size(zz,1));
str = {};
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
        t_fbest = xx(structure.tf_vars);
        x0_best = x_0;
        x0_best(~structure.imposed_initial_states) = xx(structure.x0_vars);
        xx = xx(~structure.static_vars);
        [x_best,u_best,x_b] = extract_solution(xx,structure,x_f);
        %qq = structure.uniform_in_nodes_state*t_fbest;
        %time = sort(qq(:));
        time = linspace(0,t_fbest,1000);
        [xt,ut] = eval_solution_over_time(x_best,u_best,0,t_fbest,time,structure.uniform_els,structure);
        
        % get individual states over time
        r = xt(:,1);
        vr = xt(:,2);
        theta = xt(:,3);
        vt = xt(:,4);
        
        % get individual controls over time
        u1 = ut(:,1);

        figure(1) % x over t
        box on
        set(gca,'LooseInset',get(gca,'TightInset'));
        plot(time,r,'Color',Colorset(i,:));
        hold on
        xlabel('time [TU]');
        ylabel('r [LU]');
        
        figure(2) % vx over t
        box on
        set(gca,'LooseInset',get(gca,'TightInset'));
        plot(time,vr,'Color',Colorset(i,:));
        hold on
        xlabel('time [TU]');
        ylabel('v_r [LU/TU]');
        
        figure(3) % y over t
        box on
        set(gca,'LooseInset',get(gca,'TightInset'));
        plot(time,theta,'Color',Colorset(i,:));
        hold on
        xlabel('time [TU]');
        ylabel('\theta [rad]');
        
        figure(4) % vy over t
        box on
        set(gca,'LooseInset',get(gca,'TightInset'));
        plot(time,vt,'Color',Colorset(i,:));
        hold on
        xlabel('time [TU]');
        ylabel('v_t [LU/TU]');
        
        figure(5) % q1 over t
        box on
        set(gca,'LooseInset',get(gca,'TightInset'));
        plot(time,u1,'Color',Colorset(i,:));
        hold on
        xlabel('time [TU]');
        ylabel('u [rad]');
        
        figure(6)
        box on
        set(gca,'LooseInset',get(gca,'TightInset'));
        plot(r.*cos(theta),r.*sin(theta),'Color',Colorset(i,:));
        hold on
        xlabel('x [LU]');
        ylabel('y [LU]');        
        
        thisstr = ['Solution ', num2str(i)];
        str = [str; thisstr];
        
    end
    
    %str = {'Solution 1','Solution 2','Solution 3','Solution 4','Solution 5','Solution 6','Solution 7','Solution 8','Solution 9','Solution 10'};
    
    for i = 1:6
        
        figure(i)
        legend(str)
        
    end
    
    % plot 'x' to mark the end of a spiral
    for i = 1:size(zz,1)
        
        % get data
        xx = zz(i,1:length(x_guess));
        t_fbest = xx(structure.tf_vars);
        x0_best = x_0;
        x0_best(~structure.imposed_initial_states) = xx(structure.x0_vars);
        xx = xx(~structure.static_vars);
        [x_best,u_best,x_b] = extract_solution(xx,structure,x_f);
        %qq = structure.uniform_in_nodes_state*t_fbest;
        %time = sort(qq(:));
        time = linspace(0,t_fbest,1000);
        [xt,ut] = eval_solution_over_time(x_best,u_best,0,t_fbest,time,structure.uniform_els,structure);
        
        % get individual states over time
        r = xt(:,1);
        theta = xt(:,3);

        figure(6)
        plot(r(end).*cos(theta(end)),r(end).*sin(theta(end)),'x','Color',Colorset(i,:));
        
    end
    
    
end