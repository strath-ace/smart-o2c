close all
clear all
clc

%% Equations, initial conditions and time span

f = @Simple_state_equations;
dfx = [];%@Simple_state_jacobian;
dfu = [];%@Simple_control_jacobian;

t_0 = 0;

x_0 = [2e4 0 500 0 300e3];              % [h=20km, theta=0, v=500m/s, gamma=0, m=300ton]

imposed_final_states = [1 0 0 1 0];     % mask vector with the variables on which a final condition is imposed
x_f = [4e4 0 500 0 0];                  % vector of final conditions


%% Discretisation setting

num_elems = 4;
state_order = 5;
control_order = 5;
DFET = 1;
state_distrib = 'Lobatto';
control_distrib = 'Legendre';

integr_order = 2*state_order;%2*state_order-1;

num_eqs = length(x_0);
num_controls = 2;

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

%% Optimisation

state_bounds = [0 5e4; -pi pi; 1 10000; -pi/2 pi/2; 50e3 300e3];              % h, theta, v, gamma, m
control_bounds = [0 1; -10 50];                                         % delta (throttle) from 0 to 1, alpha (angle of attack) from -10 to 50 (DEGREES, INTERNAL DYNAMICS CONVERTS CONSISTENTLY INTO RADIANS WHEN NEEDED)

[lbv,ubv] = transcribe_bounds(state_bounds,control_bounds,structure);

% initial guess
t_min = 0.5*60;%(max(state_bounds(5,:))-min(state_bounds(5,:)))/850;    % t_min mission time, considering constant max fuel flow equal to 850kg/s
t_max = 1*60;                                                  % max mission time, 1hour, maybe better estimate helps
t_guess = t_min;

lbv = [t_min;lbv];
ubv = [t_max;ubv];

% normalisation stuff

structure.scale_primal_states = state_bounds(:,2)-state_bounds(:,1);
structure.offset_primal_states = state_bounds(:,1);
structure.scale_primal_controls = control_bounds(:,2)-control_bounds(:,1);
structure.offset_primal_controls = control_bounds(:,1);
structure.scale_optimisation_vars = ubv-lbv;
structure.offset_optimisation_vars = lbv;
structure.scale_transcription_vars = structure.scale_optimisation_vars(2:end);
structure.offset_transcription_vars = structure.offset_optimisation_vars(2:end);
structure.scale_tf = structure.scale_optimisation_vars(1);
structure.offset_tf = 0;

control_guess = [0.9 5];                                        % 90% throttle and 2 degrees angle of attack
norm_control = (control_guess-structure.offset_primal_controls')./structure.scale_primal_controls'; % normalised

norm_t_guess = (t_guess-t_min)/(t_max-t_min);
norm_x0 = (x_0-structure.offset_primal_states')./structure.scale_primal_states';
norm_xf = (x_f-structure.offset_primal_states')./structure.scale_primal_states';
norm_lbv = zeros(size(lbv));
norm_lbv(1) = t_min/t_max;
norm_ubv = ones(size(ubv));

u_nodes = [norm_control(1)*ones(num_elems*(control_order+1),1) norm_control(2)*ones(num_elems*(control_order+1),1)]; % normalised

x_guess = make_first_guess(f,norm_x0,t_0,norm_t_guess,u_nodes,structure);
x_guess = [norm_t_guess;x_guess];                                      % t_final is a static optimisation parameter

%utyp = [2*t_min; utyp];
 
%options = optimset('Display','iter','GradConstr','on','MaxFunEvals',100000);
options = optimoptions(@fmincon,'Algorithm','sqp','Display','iter-detailed','MaxFunEvals',1000000,'MaxIter',1000000,'GradConstr','on','TolFun',1e-6,'Tolcon',1e-6,'TolX',1e-20,'DerivativeCheck','off');

%x_sol = fmincon(@(x) 1,x_guess,[],[],[],[],lbv,ubv,@(x) Simple_dynamics(f,structure,x,x_0,x_f,t_0,dfx,dfu),options);
[x_sol,~,~,~,~,grad] = fmincon(@(x) Simple_objectives(x,t_0,norm_xf,structure),x_guess,[],[],[],[],norm_lbv,norm_ubv,@(x) Simple_dynamics(f,structure,x,norm_x0,norm_xf,t_0,dfx,dfu),options);
x_sol = x_sol.*structure.scale_optimisation_vars+structure.offset_optimisation_vars;
t_fbest = x_sol(1);
xx = x_sol(2:end);
[x_best,u_best,x_b] = extract_solution(xx,structure,x_f);

%% Plot best solution

plot_solution_vs_time(x_best,u_best,x_0,x_b,t_0,t_fbest,structure);
% tplot = linspace(t_0,t_fbest,t_fbest*10);
% [xt,ut] = eval_solution_over_time(x_best,u_best,t_0,t_fbest,tplot,structure);
% 
% for i =1:size(xt,1)
%    
%     figure(1)
%     plot((xt(1:i,1)+3680e3).*sin(xt(1:i,2)),xt(1:i,1)+3680e3,'b');
%     hold on
%     quiver((xt(i,1)+3680e3).*sin(xt(i,2)),(xt(i,1)+3680e3),xt(i,3)*cos(xt(i,4)),xt(i,3)*sin(xt(i,4)),'Color','r','AutoScaleFactor',1e1);
%     axis([0 60e3 3680e3 3680e3+40e3])
%     drawnow
%     hold off
%     
% end