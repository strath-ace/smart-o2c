close all
clear all
clc

%% Equations, initial conditions and time span

f = @Simple_state_equations;
dfx = [];%@Simple_state_jacobian;
dfu = [];%@Simple_control_jacobian;

t_0 = 0;

x_0 = [2e4 0 500 0 300e3];              % [h=20km, theta=0, v=500m/s, gamma=0, m=300ton]

imposed_final_states = [1 0 1 1 0];     % mask vector with the variables on which a final condition is imposed
x_f = [4e4 0 500 0 0];                  % vector of final conditions


%% Discretisation setting

num_elems = 1;
state_order = 5;
control_order = 5;
DFET = 1;
state_distrib = 'Legendre';
control_distrib = 'Legendre';

integr_order = 2*(state_order);%2*state_order-1;
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

state_bounds = [0 5e4; -pi pi; 1 1000; -pi/2 pi/2; 50e3 350e3];              % h, theta, v, gamma, m
control_bounds = [0 1; -10 50];                                         % delta (throttle) from 0 to 1, alpha (angle of attack) from -10 to 50 (DEGREES, INTERNAL DYNAMICS CONVERTS CONSISTENTLY INTO RADIANS WHEN NEEDED)

[lbv,ubv] = transcribe_bounds(state_bounds,control_bounds,structure);

% initial guess
t_min = 0.5*60;%(max(state_bounds(5,:))-min(state_bounds(5,:)))/850;    % t_min mission time, considering constant max fuel flow equal to 850kg/s
t_max = 10*60;                                                  % max mission time, 1hour, maybe better estimate helps
t_guess = 4*t_min;

min_els_size = 1/5;%0.05;
minels = min_els_size+min_els_size*linspace(0,1,num_elems-1)';
maxels = flipud(1-minels);%[0.7;0.8;0.9];

lbv = [t_min;lbv];
ubv = [t_max;ubv];


%lbv = [t_min;minels;lbv];
%ubv = [t_max;maxels;ubv];

% normalisation stuff

structure.scale_primal_states = state_bounds(:,2)-state_bounds(:,1);
structure.offset_primal_states = state_bounds(:,1);
structure.scale_primal_controls = control_bounds(:,2)-control_bounds(:,1);
structure.offset_primal_controls = control_bounds(:,1);
structure.scale_optimisation_vars = ubv-lbv;
structure.offset_optimisation_vars = lbv;
structure.scale_transcription_vars = structure.scale_optimisation_vars(2+num_elems-1:end);
structure.offset_transcription_vars = structure.offset_optimisation_vars(2+num_elems-1:end);
structure.scale_tf = structure.scale_optimisation_vars(1);
structure.offset_tf = 0;

control_guess = [0.1 5];                                        % 90% throttle and 2 degrees angle of attack
norm_control = (control_guess-structure.offset_primal_controls')./structure.scale_primal_controls'; % normalised

norm_t_guess = (t_guess-t_min)/(t_max-t_min);
norm_x0 = (x_0-structure.offset_primal_states')./structure.scale_primal_states';
norm_xf = (x_f-structure.offset_primal_states')./structure.scale_primal_states';

norm_lbv = [zeros(size(structure.scale_optimisation_vars,1),1)];
norm_lbv(1) = t_min/t_max;
norm_ubv = [ones(size(structure.scale_optimisation_vars,1),1)];


% norm_lbv = [t_min/t_max;minels;zeros(size(structure.scale_transcription_vars,1),1)];
% norm_ubv = [1;maxels;ones(size(structure.scale_transcription_vars,1),1)];

u_nodes = [norm_control(1)*ones(num_elems*(control_order+1),1) norm_control(2)*ones(num_elems*(control_order+1),1)]; % normalised

x_guess = make_first_guess(f,norm_x0,t_0,norm_t_guess,u_nodes,structure);
x_guess = [norm_t_guess;x_guess];                                      % t_final is a static optimisation parameter
%x_guess = [norm_t_guess;structure.uniform_els(2:end,1);x_guess];                                      % t_final is a static optimisation parameter

options = optimoptions(@fmincon,'Algorithm','sqp','Display','iter-detailed','MaxFunEvals',1000000,'MaxIter',1000000,'GradConstr','off','TolFun',1e-6,'Tolcon',1e-6,'TolX',1e-20,'DerivativeCheck','off');

b = -min_els_size*ones((2+num_elems-2)*(num_elems>1),(num_elems>1));
A = zeros((2+num_elems-2)*(num_elems>1),length(x_guess));                                 % Inequality Constraints between elements (size)

if size(A,1)>0
    
    b(end) = 1-min_els_size;
    A(1,2) = -1;
    
    for i = 2:num_elems-1
        
        A(i,i:i+1) = [1 -1];
        
    end
    
    A(end,num_elems)= 1;
    
end

[x_sol,~,~,~,~,grad] = fmincon(@(x) Simple_objectives(x,t_0,norm_xf,structure),x_guess,[],[],[],[],norm_lbv,norm_ubv,@(x) Simple_dynamics(f,structure,x,norm_x0,norm_xf,t_0,dfx,dfu),options);
%[x_sol,~,~,~,~,grad] = fmincon(@(x) Simple_objectives(x,t_0,norm_xf,structure),x_guess,A,b,[],[],norm_lbv,norm_ubv,@(x) Simple_dynamics(f,structure,x,norm_x0,norm_xf,t_0,dfx,dfu),options);
x_sol = x_sol.*structure.scale_optimisation_vars+structure.offset_optimisation_vars;
t_fbest = x_sol(1);
els_best = structure.uniform_els;
%els_best = x_sol(2:2+structure.num_elems-2);
%els_best = [[0; els_best] [els_best;1]];
%xx = x_sol(2+structure.num_elems-1:end);
xx = x_sol(2:end);
[x_best,u_best,x_b] = extract_solution(xx,structure,x_f);

%% Plot best solution

plot_solution_vs_time(x_best,u_best,x_0,x_b,t_0,t_fbest,els_best,structure);

%plot_solution_vs_time(x_best,u_best,x_0,x_b,t_0,t_fbest,structure.uniform_els,structure);

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