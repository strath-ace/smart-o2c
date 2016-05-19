close all
clear all
clc

%% Equations, initial conditions and time span

f = @state_equation;
dfx = @state_jacobian;
dfu = @control_jacobian;

t_0 = 0;

x_0 = [2e4 0 500 0 300e3];              % [h=20km, theta=0, v=500m/s, gamma=0, m=300ton]

imposed_final_states = [1 0 0 0 0];     % mask vector with the variables on which a final condition is imposed
x_f = [4e4 0 500 0 0];                  % vector of final conditions


%% Discretisation settings

num_elems = 50;
state_order = 0;
control_order = 0;
DFET = 1;
state_distrib = 'Lobatto'; % or Cheby
control_distrib = 'Legendre';

%state_nodes_distrib = 'Legendre';
%control_nodes_distrib = 'Legendre';

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

state_bounds = [0 5e4; 0 2*pi; 1 1000; -pi/2 pi/2; 50e3 300e3];              % h, theta, v, gamma, m
control_bounds = [0 1; -10 50];                                         % delta (throttle) from 0 to 1, alpha (angle of attack) from -10 to 50 (DEGREES, INTERNAL DYNAMICS CONVERTS CONSISTENTLY INTO RADIANS WHEN NEEDED)

[lbv,ubv] = transcribe_bounds(state_bounds,control_bounds,structure);

%typical values (experiment)
%state_typical = [0 2e4; 0 2*pi; 0 500; 0 2*pi; 50e3 300e3];
%control_typical = control_bounds;

%[ltyp,utyp] = transcribe_bounds(state_typical,control_typical,structure);

% initial guess
u_nodes = [ones(num_controls*num_elems*(control_order+1),1) 2*ones(num_controls*num_elems*(control_order+1),1)]; % half throttle and 10 degrees angle of attack

t_min = 0.5*60;%(max(state_bounds(5,:))-min(state_bounds(5,:)))/850;    % t_min mission time, considering constant max fuel flow equal to 850kg/s
t_max = 10*60;                                                  % max mission time, 1hour, maybe better estimate helps
t_guess = 2*t_min;

x_guess = make_first_guess(f,x_0,t_0,t_guess,u_nodes,structure);
x_guess = [t_guess;x_guess];                                      % t_final is a static optimisation parameter

lbv = [t_min;lbv];
ubv = [t_max;ubv];

%utyp = [2*t_min; utyp];
 
%options = optimset('Display','iter','GradConstr','on','MaxFunEvals',100000);
options = optimoptions(@fmincon,'Algorithm','sqp','Display','iter-detailed','MaxFunEvals',1000000,'MaxIter',1000000,'GradConstr','on','TolFun',1e-6,'TolX',1e-15,'DerivativeCheck','off');

%x_sol = fmincon(@(x) 1,x_guess,[],[],[],[],lbv,ubv,@(x) dynamics(f,structure,x,x_0,x_f,t_0,dfx,dfu),options);
x_sol = fmincon(@(x) objectives(x,t_0,x_f,structure),x_guess,[],[],[],[],lbv,ubv,@(x) dynamics(f,structure,x,x_0,x_f,t_0,dfx,dfu),options);
t_fbest = x_sol(1);
xx = x_sol(2:end);
[x_best,u_best,x_b] = extract_solution(xx,structure,x_f);

%% Plot best solution

plot_solution_vs_time(x_best,u_best,x_0,x_b,t_0,t_fbest,structure);