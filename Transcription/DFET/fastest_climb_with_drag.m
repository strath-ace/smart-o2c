close all
clear all
clc

maxtrials = 1;

%% Equations, initial conditions and time span

% formulation with "internal wrapping"

f = @state_equation;
dfx = @state_jacobian;
dfu = @control_jacobian;

%f = @(x,u,t) [x(2); a*cos((u(1)))-0.5*rho0*exp(-x(3))*S*Cd*(x(2).^2+x(4).^2).^(0.5)*x(2); x(4); a*sin((u(1)))-g-0.5*rho0*exp(-x(3))*S*Cd*(x(2).^2+x(4).^2).^(0.5)*x(4)];
%dfx = @(x,u,t) [0 1 0 0; 0 -0.5*rho0*exp(-x(3))*S*Cd*(x(2).^2/(x(2).^2+x(4).^2).^0.5+(x(2).^2+x(4).^2).^0.5).*((x(2).^2+x(4).^2)>0) 0.5*rho0*exp(-x(3))*S*Cd*(x(2).^2+x(4).^2).^(0.5)*x(2) -0.5*rho0*exp(-x(3))*S*Cd*x(2)*x(4)/(x(2).^2+x(4).^2).^0.5*((x(2).^2+x(4).^2)>0); 0 0 0 1; 0 -0.5*rho0*exp(-x(3))*S*Cd*x(2)*x(4)/(x(2).^2+x(4).^2).^0.5*((x(2).^2+x(4).^2)>0) 0.5*rho0*exp(-x(3))*S*Cd*(x(2).^2+x(4).^2).^(0.5)*x(4) -0.5*rho0*exp(-x(3))*S*Cd*(x(4).^2/(x(2).^2+x(4).^2).^0.5+(x(2).^2+x(4).^2).^0.5)*((x(2).^2+x(4).^2)>0)];
%dfu = @(x,u,t) [0; -a*sin(u(1)); 0; a*cos(u(1))];
%dfx = [];
%dfu = [];

t_0 = 0;

x_0 = [0 0 0 0 1];

imposed_final_states = [0 0 1 1 0];   % mask vector with the variables on which a final condition is imposed
x_f = [0 0.1 10 0 0];                % vector of final conditions (size???)


%% Discretisation settings

num_elems = 100;
state_order = 4;
control_order = 0;
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

%% Transcribe constraints

structure = prepare_transcription(num_eqs,num_controls,num_elems,state_order,control_order,integr_order,DFET,state_distrib,control_distrib);

structure = impose_final_conditions(structure,imposed_final_states);

%% Optimisation loop

best_cost = inf(maxtrials,1);
costs = nan(maxtrials,1);

state_bounds = [-inf inf;-inf inf;-inf inf;-inf inf; 0 inf];
control_bounds = [-pi/2 pi/2];

[lbv,ubv] = transcribe_bounds(state_bounds,control_bounds,structure);

u_nodes = pi/4*ones(num_controls*num_elems*(control_order+1),1);
%u_nodes = linspace(0.6,0.2,length(u_nodes))';
x_guess = make_first_guess(f,x_0,t_0,250,u_nodes,structure);
x_guess = [120;x_guess];

lbv = [10;lbv];
ubv = [250;ubv];
 
%options = optimset('Display','iter','GradConstr','on','MaxFunEvals',100000);
options = optimset('Display','iter-detailed','MaxFunEvals',1000000,'MaxIter',1000000,'GradConstr','on','TolFun',1e-8,'TolX',1e-12);

x_sol = fmincon(@(x) objectives(x,t_0,x_f,structure),x_guess,[],[],[],[],lbv,ubv,@(x) dynamics(f,structure,x,x_0,x_f,t_0,dfx,dfu),options);
t_fbest = x_sol(1);
xx = x_sol(2:end);
% t_fbest = 250;%x_sol(1);
% xx = x_sol;%(2:end);
[x_best,u_best,x_b] = extract_solution(xx,structure,x_f);
%u_best = wrapToPi(u_best);

% maxtrials = 250;
% times = linspace(250,108,maxtrials);
% best_cost = inf(maxtrials,1);
% 
%  for i=1:maxtrials
%     
%     t_f = times(i);
%     
%     %t_f = 250+40*rand();
%     
%     g = @(x,u,t) [t_f 0];
%     weights = [1 0];
%     
%     % bogus initial control parameters
%     
%     u_nodes = pi/2*ones(num_controls*num_elems*(control_order+1),1);
%     
%     % Initial guess
%     
%     x_guess = make_first_guess(f,x_0,t_0,t_f,u_nodes,structure);
%     
%     % Solution of Under-constrained Nonlinear System
%     tic
%     [x_sol,resids] = solve_constraints(f,x_guess,x_0,x_f,t_0,t_f,structure,1e-3,50000,dfx,dfu,lbv,ubv);
%     toc
%     % Extract solutions (separe state/controls)
%     
%     [x,u,x_b] = extract_solution(x_sol,structure,x_f);
%     
%     % Evaluation of cost function
%     costs(i) = eval_cost_functions(g,weights,x,u,x_b,[t_0 t_f],structure);
%     
%     costs(i) = costs(i)*(resids(end)<1e-3)+250*(resids(end)>=1e-3);
%     
%     if costs(i)<best_cost(i)
%         
%         best_cost(i:end) = costs(i);
%         x_best = x;
%         u_best = u;
%         best_str = structure;
%         t_fbest = t_f;
%         
%     end
%     
%     plot(1:i,best_cost(1:i),'b*',1:i,costs(1:i),'ro')
%     drawnow
%     
% end

%% Plot best solution

% wrap controls between -pi/2 and pi/2
%u_best = atan(sin(u_best)./cos(u_best));

% other way to wrap controls between -pi/2 and pi/2, maybe less problem
% dependent
 
% while max(max(u_best))>pi/2
%     
%     u_best = u_best.*(u_best<=pi/2) + (u_best-2*pi).*(u_best>pi/2);
%     
% end
% 
% while min(min(u_best))<-pi/2
%     
%     u_best = u_best.*(u_best>=-pi/2) + (u_best+2*pi).*(u_best<-pi/2);
%     
% end


plot_solution_vs_time(x_best,u_best,x_0,x_b,t_0,t_fbest,structure);
tplot = linspace(t_0,t_fbest,100);
 %tplot = structure.in_nodes_state';
 %tplot = tplot(:);
 %tplot = tplot'*t_fbest;
[xt,ut] = eval_solution_over_time(x_best,u_best,t_0,t_fbest,tplot,structure);
figure()
plot(tplot,cos(ut),'b.')
figure()
plot(tplot,sin(ut),'b.')

% if DFET==1
%    
%     errs = [xt(1,:)-x_0; xt(2:test_order-1:end,:)-xt(1:test_order-1:end-1,:); xt(end,:)-x_b]
%     
% end