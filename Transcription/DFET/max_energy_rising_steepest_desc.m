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
t_f = 10;

x_0 = [1.1 0 0 1/(1.1)^0.5];

imposed_final_states = [0 0 0 0];   % mask vector with the variables on which a final condition is imposed
x_f = [0 0 0 0];                % vector of final conditions (size???)

% cost functions and weights
g = @(x,u,t) [-(1/2*(x(2)^2+x(4)^2))-1/x(1) 0];
dgxf = @(x,u,t) [-1/x(1)^2 -x(2) 0 -x(4)];
dguf = @(x,u,t) [0];
dgxi = @(x,u,t) [0 0 0 0];
dgui = @(x,u,t) [0];
weights = [1 0];

%% Discretisation settings

num_elems = 10;
state_order = 1;
control_order = 1;
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

u_nodes = 4*pi/9*ones(num_controls*num_elems*(control_order+1),1);

% Initial guess

x_guess = make_first_guess(f,x_0,t_0,t_f,u_nodes,structure);

% Solution of Under-constrained Nonlinear System

[x_sol,resids] = solve_constraints(f,x_guess,x_0,x_f,t_0,t_f,structure,1e-6,5000,dfx,dfu,lbv,ubv);

% Extract solutions (separe state/controls)

[x,u,x_b] = extract_solution(x_sol,structure,x_f);

% Evaluation of cost function
[vals,gradu] = eval_cost_functions2(g,weights,x,u,x_b,[t_0 t_f],structure,1,dfu,dgxf,dguf,dgxi,dgui);

x_best = x;
u_best = u;
x_b_best = x_b;
g_best = vals;
gradu_best = gradu;
x_sol_best = x_sol;
figure(1)
plot_solution_vs_time(x_best,u_best,x_0,x_b_best,t_0,t_f,structure);

fprintf('Iteration \t Objective_fun_value \t Step_size \t Optimality \t Residual \n');
fprintf('%d \t\t %f \t\t %s \t\t %f \t %f \n',0,vals,'*',max(abs(gradu)),resids(end));

maxit = 200*control_order*num_elems*num_controls;
it = 1;

maxtrials = 200*(control_order+1)*num_controls*num_elems;
s = 0.1;

max_fails = 10;
fails = 0;
for i=1:maxtrials
        
    u_nodes = u_best(:)-s*gradu_best';
    
    if ~isempty(control_bounds)
        
        [s,gradu] = alpha_clip2(u_best(:),s,gradu',repmat(control_bounds(1),length(u_nodes),1),repmat(control_bounds(2),length(u_nodes),1));
    
        s = s*0.99; % should avoid roundoff problems
        if any((u_best(:)+s*gradu)>repmat(control_bounds(2),length(u_nodes),1)) || any((u_best(:)+s*gradu)<repmat(control_bounds(1),length(u_nodes),1))
            
            error('failed to constrain bounds!!!');
            
        end
        
    end  
    
    x_guess = update_guess(u_nodes,x_best,structure);
    
    % Solution of Under-constrained Nonlinear System
    
    [x_sol,resids] = solve_constraints(f,x_guess,x_0,x_f,t_0,t_f,structure,1e-6,5000,dfx,dfu,lbv,ubv);
    
    % Extract solutions (separe state/controls)
    
    [x,u,x_b] = extract_solution(x_sol,structure,x_f);
    
    if resids(end)<1e-6
        
        % Extract solutions (separe state/controls)
        [x,u,x_b] = extract_solution(x_sol,structure,x_f);
        
        % Eval cost function and their gradients
        [vals,gradu] = eval_cost_functions2(g,weights,x,u,x_b,[t_0 t_f],structure,1,dfu,dgxf,dguf,dgxi,dgui);
        
    else
        
        vals = nan;
        
    end
    
    if vals<g_best
        
        u_best = u;
        g_best = vals;
        x_best = x;
        x_sol_best = x_sol;
        gradu_best = gradu;
        
        % Plot solution
        subplot(1,1,1)
        cla reset
        subplot(2,1,2)
        cla reset
        plot_solution_vs_time(x,u,x_0,x_b,t_0,t_f,structure);
        drawnow
        fprintf('%d \t\t %f \t\t %f \t %f \t %f \n',it,vals,s,max(abs(gradu)),resids(end));
        s=1;
        
    else
        
        fprintf('%d \t\t %f \t\t %f \t %f \t (%f)* \n',it,g_best,s,max(abs(gradu_best)),resids(end));

        s = s*0.5;
        
    end
    
%     if fails==max_fails
%        
%         s=s*0.5;
%         fails = 0;
%         
%     end
    
    
    it = it+1;
    
end

% Plot best trajectory

tplot = linspace(t_0,t_f,100);
[xt,ut] = eval_solution_over_time(x_best,u_best,t_0,t_f,tplot,structure);

figure(2)
plot(xt(:,1).*cos(xt(:,3)),xt(:,1).*sin(xt(:,3)),'b')
axis equal
hold on
plot(0.5*cos([0:0.01:2*pi]),0.5*sin([0:0.01:2*pi]),'k')
