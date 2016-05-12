close all
clear all
clc

maxtrials = 1;
a = 4e-3;

%% Equations, initial conditions and time span

% formulation with "internal wrapping"
%f = @(x,u,t) [x(2); a*cos(wrapToPi(u(1))); x(4); -0.0016+a*sin(wrapToPi(u(1)))];

f = @(x,u,t) [x(2); a*cos((u(1))); x(4); -0.0016+a*sin((u(1)))];

dfx = @(x,u,t) [0 1 0 0; 0 0 0 0; 0 0 0 1; 0 0 0 0];
dfu = @(x,u,t) [0; -a*sin(u(1)); 0; a*cos(u(1))];
%dfx = [];
%dfu = [];

t_0 = 0;

x_0 = [0 0 0 0];

imposed_final_states = [0 1 1 1];   % mask vector with the variables on which a final condition is imposed
x_f = [0 0.1 10 0];                % vector of final conditions (size???)

% cost functions and weights
t_f=120;
g = @(x,u,t) [t_f 0];
weights = [1 0];

%% Discretisation settings

num_elems = 4;
state_order = 8;
control_order = 8;
DFET = 1;
state_distrib = 'Lobatto'; % or Cheby
control_distrib = 'Legendre';

%state_nodes_distrib = 'Legendre';
%control_nodes_distrib = 'Legendre';

integr_order = 2*state_order;%2*state_order-1;

num_eqs = length(x_0);
num_controls = 1;

umax = 1;

%% Make checks

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

%% Transcribe problem and box constraints

structure = prepare_transcription(num_eqs,num_controls,num_elems,state_order,control_order,integr_order,DFET,state_distrib,control_distrib);

structure = impose_final_conditions(structure,imposed_final_states);

state_bounds = [-inf inf;-inf inf;-inf inf;-inf inf];
control_bounds = [-pi/2 pi/2];

[lbv,ubv] = transcribe_bounds(state_bounds,control_bounds,structure);

%% Solution

% bogus initial control parameters

u_nodes = pi/2*ones(num_controls*num_elems*(control_order+1),1);

% Initial guess
x_guess = make_first_guess(f,x_0,t_0,t_f,u_nodes,structure);

% Solution of Under-constrained Nonlinear System
tolconv = 1e-6;
maxit_solv = 5000;
[x_sol,resids] = solve_constraints(f,x_guess,x_0,x_f,t_0,t_f,structure,tolconv,maxit_solv,dfx,dfu,lbv,ubv);

% Extract solutions (separe state/controls)

[x,u,x_b] = extract_solution(x_sol,structure,x_f);

[vals] = eval_cost_functions(g,weights,x,u,x_b,[t_0 t_f],structure,0,[],[],[],[]);

plot_solution_vs_time(x,u,x_0,x_b,t_0,t_f,structure);
drawnow

fprintf('Iteration \t Objective_fun_value \t Step_size \t Residual \n');
fprintf('%d \t\t %e \t\t\t %s \t\t\t %e \n',0,vals,'*',resids(end));

maxit = 200*control_order*num_elems*num_controls;
minstep = 1e-6;
it = 1;

u_best = u;
g_best = vals;
x_best =  x;
t_f_best = t_f;

s = 1;
%x_guess = x_sol;

while it<=maxit && s>minstep
    
    t_f = t_f_best-s*rand();%+2*rand();
    
    if t_f>150
       
        t_f = 150;
        
    end
    
    if t_f<100
       
        t_f = 100;
        
    end
    
    g = @(x,u,t) [t_f 0];
    
    perturb = -1+2*rand(length(u_nodes),1);
    
    u_nodes = u_best(:)+s*perturb;

    u_nodes(u_nodes>pi/2) = pi/2;
    u_nodes(u_nodes<-pi/2) = -pi/2;

    x_guess = update_guess(u_nodes,x_best,structure);
    
    % Solution of Under-constrained Nonlinear System
    [x_sol,resids,tt] = solve_constraints(f,x_guess,x_0,x_f,t_0,t_f,structure,tolconv,maxit_solv,dfx,dfu,lbv,ubv);
    
    if resids(end)<1e-6
        
        % Extract solutions (separe state/controls)
        [x,u,x_b] = extract_solution(x_sol,structure,x_f);
        
        % Eval cost function and their gradients
        [vals] = eval_cost_functions(g,weights,x,u,x_b,[t_0 t_f],structure,0,[],[],[],[]);
        
    else
        
        vals = nan;
        
    end
    
    if vals<g_best
        
        u_best = u;
        g_best = vals;
        x_best = x;
        t_f_best = t_f;
        
        % Plot solution
        figure(1)
        subplot(1,1,1)
        cla reset
        subplot(2,1,2)
        cla reset
        plot_solution_vs_time(x,u,x_0,x_b,t_0,t_f,structure);
        fprintf('%d \t\t %e \t\t\t %e \t %e \n',it,vals,s,resids(end));
        %s=1;
        
    else
        
        fprintf('%d \t\t %e \t\t\t %e \t (%e)* \n',it,g_best,s,resids(end));
        s=s*0.5;
        
    end
    
    figure(2)
    cla reset
    loglog(1:length(resids),resids,1:length(resids),tolconv*ones(length(resids),1));
    hold on
    semilogx(1:length(resids),tt,'g');
    drawnow
    it = it+1;
    
end
