close all
clear all
clc

%% Equations, initial conditions and time span


f = @(x,u,t) [x(2); 1.4*x(2)-0.14*x(2).^3-x(1)+4*u(1)];
dfx = @(x,u,t) [0 1; -1 1.4-0.14*3*x(2)^2];
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


%% Discretisation settings

num_elems = 1;
state_order = 20;
control_order = 1;
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
tic
[x_sol,resids] = solve_constraints(f,x_guess,x_0,x_f,t_0,t_f,structure,1e-6,500,dfx,dfu,[],[]);
toc
% Extract solutions (separe state/controls)

[x,u,x_b] = extract_solution(x_sol,structure,x_f);

%% Eval cost function
g = @(x,u,t) [0 u(1)^2; u(1)^2 0];

% order is: derivatives of first function (all variables), THEN derivatives
% of integral of second function (all variables)
dgxf = @(x,u,t) [0 0; 0 0];
dguf = @(x,u,t) [0; 2*u(1)];
dgxi = @(x,u,t) [0 0; 0 0];
dgui = @(x,u,t) [2*u(1); 0];

weights = [0 1; 1 0];
[vals,gradx,gradu] = eval_cost_functions(g,weights,x,u,x_b,[t_0 t_f],structure,1,dgxf,dguf,dgxi,dgui);

%% Plot best solution

% plot_solution_vs_time(x,u,x_0,x_b,t_0,t_f,structure);
% tplot = linspace(t_0,t_f,100);
% [xt,ut] = eval_solution_over_time(x,u,t_0,t_f,tplot,structure);
% 
% if DFET==1
%    
%     errs = [xt(1,:)-x_0; xt(end,:)-x_b]'
%     
% end