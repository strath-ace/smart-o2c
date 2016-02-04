close all
clear all
clc

maxtrials = 1;
a = 1e-2;

%% Equations, initial conditions and time span

f = @(x,u,t) [x(2); x(4)^2/x(1)-1/x(1)^2+a*cos(u(1)); x(4)/x(1); -x(2)*x(4)/x(1)+a*sin(u(1))];
dfx = @(x,u,t) [0 1 0 0; -x(4)^2/x(1)^2+2/(x(1)^3) 0 0 2*x(4)/x(1); -x(4)/(x(1)^2) 0 0 1/x(1); x(2)*x(4)/(x(1)^2) -x(4)/x(1) 0 -x(2)/x(1)];
dfu = @(x,u,t) [0; -a*sin(u(1)); 0; a*cos(u(1))];

t_0 = 0;
t_f = 20;

x_0 = [1.1 0 0 1/(1.1)^0.5];

imposed_final_states = [0 0 0 0];   % mask vector with the variables on which a final condition is imposed
x_f = [0 0 0 0];                % vector of final conditions (size???)

% cost functions and weights
g = @(x,u,t) [-(1/2*(x(2)^2+x(4)^2)-1/x(1)) 0];
weights = [1 0];

%% Discretisation settings

num_elems = 30;
state_order = 3;
control_order = 3;
DFET = 1;
state_distrib = 'Lobatto'; % or Cheby
control_distrib = 'Legendre';

integr_order = 2*state_order;%2*state_order-1;

num_eqs = length(x_0);
num_controls = 1;

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

%% Transcribe constraints

structure = prepare_transcription(num_eqs,num_controls,num_elems,state_order,control_order,integr_order,DFET,state_distrib,control_distrib);

structure = impose_final_conditions(structure,imposed_final_states);

%% Optimisation loop

best_cost = inf(maxtrials,1);
costs = nan(maxtrials,1);

%state_bounds = [0 inf;-inf inf; -inf inf;-inf inf];
%control_bounds = [-pi pi];

state_bounds = [0 inf;-inf inf; -inf inf;-inf inf];
control_bounds = [-pi pi];

[lbv,ubv] = transcribe_bounds(state_bounds,control_bounds,structure);

% bogus initial control parameters

u_nodes = pi/2*ones(num_controls*num_elems*(control_order+1),1);

% Initial guess

x_guess = make_first_guess(f,x_0,t_0,t_f,u_nodes,structure);

% Solution of Under-constrained Nonlinear System

[x_sol,resids] = solve_constraints(f,x_guess,x_0,x_f,t_0,t_f,structure,1e-6,500,dfx,dfu,lbv,ubv);

% Extract solutions (separe state/controls)

[x,u,x_b] = extract_solution(x_sol,structure,x_f);

% Evaluation of cost function
vals = eval_cost_functions(g,weights,x,u,x_b,[t_0 t_f],structure,0,[],[],[],[]);

x_best = x;
u_best = u;
x_b_best = x_b;
g_best = vals;
x_sol_best = x_sol;
plot_solution_vs_time(x_best,u_best,x_0,x_b_best,t_0,t_f,structure);

fprintf('Iteration \t Objective_fun_value \t Step_size \t Residual \n');
fprintf('%d \t\t %e \t\t %s \t\t %e \n',0,vals,'*',resids(end));

maxtrials = 200*(control_order+1)*num_controls*num_elems;
s = 1;
it = 0;
minstep = 1e-6;

max_fails = 10;
fails = 0;

while it<maxtrials && s>minstep
    
    perturb = -1+2*rand(length(u_nodes),1);
    
    u_nodes = u_best(:)+s*perturb;
    
    u_nodes(u_nodes>pi/2) = pi/2;
    u_nodes(u_nodes<-pi/2) = -pi/2;
    
    x_guess = update_guess(u_nodes,x_best,structure);
    
    % Solution of Under-constrained Nonlinear System
    
    [x_sol,resids] = solve_constraints(f,x_guess,x_0,x_f,t_0,t_f,structure,1e-6,5000,dfx,dfu,lbv,ubv);
    
    % Extract solutions (separe state/controls)
    
    [x,u,x_b] = extract_solution(x_sol,structure,x_f);
    
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
        x_sol_best = x_sol;
        
        % Plot solution
        subplot(1,1,1)
        cla reset
        subplot(2,1,2)
        cla reset
        plot_solution_vs_time(x,u,x_0,x_b,t_0,t_f,structure);
        drawnow
        fprintf('%d \t\t %e \t\t %e \t %e \n',it,vals,s,resids(end));
        %s=1;
        fails = 0;
        
    else
        
        fprintf('%d \t\t %e \t\t %e \t (%e)* \n',it,g_best,s,resids(end));
        fails = fails+1;
        
    end
    
    if fails>5
       
        s = s*0.5;
        fails = 0;
        
    end
    
    it = it+1;
    
end

% Plot best trajectory

tplot = linspace(t_0,t_f,100);
[xt,ut] = eval_solution_over_time(x_best,u_best,t_0,t_f,tplot,structure);

figure()
plot(xt(:,1).*cos(xt(:,3)),xt(:,1).*sin(xt(:,3)),'b')
axis equal
hold on
plot(0.5*cos(tplot),0.5*sin(tplot),'k')
