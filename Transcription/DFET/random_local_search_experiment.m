close all
clear all
clc

%% Equations, initial conditions and time span


f = @(x,u,t) [x(2); 1.4*x(2)-0.14*x(2)^3-x(1)+4*u(1)];
dfx = @(x,u,t) [0 1; -1 1.4-3*0.14*x(2)^2];
dfu = @(x,u,t) [0; 4];
%dfx = @(x,u,t) [0 1 0 0; 0 0 0 0; 0 0 0 1; 0 0 0 0];
%dfu = @(x,u,t) [0; -a*sin(u(1)); 0; a*cos(u(1))];
%dfx = [];
%dfu = [];

t_0 = 0;
t_f = 2.5;

x_0 = [-5 -5];

imposed_final_states = [1 1];   % mask vector with the variables on which a final condition is imposed
x_f = [1.5 2];                % vector of final conditions (size???)

%% Definition of cost function and its analytical derivatives
g = @(x,u,t) [0 x(1)^2+u(1)^2];

% order is: derivatives of first function (all variables), THEN derivatives
% of integral of second function (all variables)
dgxf = @(x,u,t) [0 0];
dguf = @(x,u,t) [0];
dgxi = @(x,u,t) [2*x(1) 0];
dgui = @(x,u,t) [2*u(1)];

weights = [0 1];

%% Discretisation settings

num_elems = 8;
state_order = 4;
control_order = 2;
DFET = 1;
state_distrib = 'Lobatto'; % or Cheby
control_distrib = 'Legendre';

%state_nodes_distrib = 'Legendre';
%control_nodes_distrib = 'Legendre';

integr_order = max(2*state_order,1);%2*state_order-1;

num_eqs = length(x_0);
num_controls = 1;

%% Check equations/unknowns balance

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

%% Transcribe equations and constraints

structure = prepare_transcription(num_eqs,num_controls,num_elems,state_order,control_order,integr_order,DFET,state_distrib,control_distrib);

structure = impose_final_conditions(structure,imposed_final_states);

%% Solution

% bogus initial control parameters

u_nodes = 0*ones(num_controls*num_elems*(control_order+1),1);

% Initial guess

x_guess = make_first_guess(f,x_0,t_0,t_f,u_nodes,structure);

% Solution of Under-constrained Nonlinear System

[x_sol,resids] = solve_constraints(f,x_guess,x_0,x_f,t_0,t_f,structure,1e-6,500,dfx,dfu,[],[]);

% Extract solutions (separe state/controls)

[x,u,x_b] = extract_solution(x_sol,structure,x_f);

[vals] = eval_cost_functions(g,weights,x,u,x_b,[t_0 t_f],structure,1,dgxf,dguf,dgxi,dgui);

plot_solution_vs_time(x,u,x_0,x_b,t_0,t_f,structure);
drawnow

fprintf('Iteration \t Objective_fun_value\n');
fprintf('%d \t\t %f\n',0,vals);

maxit = 200*control_order*num_elems*num_controls;
it = 1;

u_best = u;
g_best = vals;
x_best =  x;

s = 1;

while it<=maxit
    
    perturb = -1+2*rand(length(u_nodes),1);
    
    u_nodes = u_best(:)+s*perturb;
    x_guess = update_guess(u_nodes,x_best,structure);
    
    % Solution of Under-constrained Nonlinear System
    [x_sol,resids] = solve_constraints(f,x_guess,x_0,x_f,t_0,t_f,structure,1e-6,500,dfx,dfu,[],[]);
    
    % Extract solutions (separe state/controls)
    [x,u,x_b] = extract_solution(x_sol,structure,x_f);
    
    % Eval cost function and their gradients
    [vals] = eval_cost_functions(g,weights,x,u,x_b,[t_0 t_f],structure,1,dgxf,dguf,dgxi,dgui);
    
    if vals<g_best
        
        u_best = u;
        g_best = vals;
        x_best = x;
        
        % Plot solution
        subplot(1,1,1)
        cla reset
        subplot(2,1,2)
        cla reset
        plot_solution_vs_time(x,u,x_0,x_b,t_0,t_f,structure);
        drawnow
        fprintf('%d \t\t %f \n',it,vals);   
        s=1;
        
    else
        
        fprintf('%d \t\t %f \n',it,g_best);   
        s=s*0.5;
    end
    
    
    it = it+1;
    
end

