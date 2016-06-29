close all
clear all
clc

%% DFET PROBLEM TRANSCRIPTION

% Add correct directory to path!!!

addpath(genpath('MO_control\Ascent\AAS'));

% Dynamics

f = @Simple_state_equations;
dfx = [];%@Simple_state_jacobian;
dfu = [];%@Simple_control_jacobian;
dft = [];
smooth_scal_const = @problem_specific_smooth_scal_constraints; 

% Objective functions and Bolza's problem weights

g = @(x,u,t) [-x(5) 0; t 0];
weights = [1 0; 1 0];

t_0 = 0;

x_0 = [2e4 0 500 0 300e3];              % [h=20km, theta=0, v=500m/s, gamma=0, m=300ton]

imposed_final_states = [1 0 1 1 0];     % mask vector with the variables on which a final condition is imposed
x_f = [4e4 0 500 0 0];                  % vector of final conditions

%% Discretisation settings

num_elems = 1;
state_order = 5;
control_order = 5;
DFET = 1;
state_distrib = 'Legendre'; % or Cheby
control_distrib = 'Legendre';

integr_order = 2*state_order;%2*state_order-1;

num_eqs = length(x_0);
num_controls = 2;

% Make checks

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

% Transcribe constraints

structure = prepare_transcription(num_eqs,num_controls,num_elems,state_order,control_order,integr_order,DFET,state_distrib,control_distrib);

structure = impose_final_conditions(structure,imposed_final_states);

% include function handles, derivatives, objective functions and weights
% for Bolza's problem in the structure. Will have to do it with a proper
% function

structure.f = f;
structure.dfx = [];
structure.dfu = [];
structure.dft = [];
structure.g = g;
structure.weights = weights;

state_bounds = [0 5e4; -pi pi; 100 1000; -pi/2 pi/2; 50e3 350e3];              % h, theta, v, gamma, m
control_bounds = [0 1; -10 50];                                         % delta (throttle) from 0 to 1, alpha (angle of attack) from -10 to 50 (DEGREES, INTERNAL DYNAMICS CONVERTS CONSISTENTLY INTO RADIANS WHEN NEEDED)

[vlb,vub,state_vars,control_vars] = transcribe_bounds(state_bounds,control_bounds,structure);

t_min = 0.5*60; 
t_max = 10*60; 

vlb = [t_min;vlb]';
vub = [t_max;vub]';

% normalisation stuff

structure.scale_primal_states = (state_bounds(:,2)-state_bounds(:,1))';
structure.offset_primal_states = state_bounds(:,1)';
structure.scale_primal_controls = (control_bounds(:,2)-control_bounds(:,1))';
structure.offset_primal_controls = control_bounds(:,1)';
structure.scale_optimisation_vars = vub-vlb;
structure.offset_optimisation_vars = vlb;
structure.scale_transcription_vars = structure.scale_optimisation_vars(2:end);
structure.offset_transcription_vars = structure.offset_optimisation_vars(2:end);
structure.scale_tf = structure.scale_optimisation_vars(1);
structure.offset_tf = t_min;

structure.scale_objectives = structure.scale_primal_states(5);%[structure.scale_primal_states(5) structure.scale_tf];
structure.offset_objectives = structure.offset_primal_states(5);%[structure.offset_primal_states(5) 0];

norm_x_0 = (x_0-structure.offset_primal_states)./structure.scale_primal_states;
norm_x_f = (x_f-structure.offset_primal_states)./structure.scale_primal_states;

structure.norm_vlb = zeros(1,length(structure.scale_optimisation_vars));
structure.norm_vlb(1) = t_min/t_max;
structure.norm_vub = ones(1,length(structure.scale_optimisation_vars));

%% MACS PARAMETERS

opt.maxnfeval=10000;                                                       % maximum number of f evals 
opt.popsize=10;                                                             % popsize (for each archive)
opt.rhoini=1;                                                               % initial span of each local hypercube (1=full domain)
opt.F=0.9;                                                                    % F, the parameter for Differential Evolution
opt.CR=0.9;                                                                   % CR, crossover probability
opt.p_social=1;                                                           % popratio
opt.max_arch=10;                                                            % archive size
opt.coord_ratio=1;%/(sum(~isinf(vlb)));                                               
opt.contr_ratio=0.5;                                                        % contraction ratio
opt.draw_flag=0;                                                            % draw flag
opt.cp=0;                                                                   % constraints yes/no 
opt.MBHflag=1;                                                              % number of MBH steps
opt.cpat=0;                                                                 % pattern to DE
opt.explore_DE_strategy = 'rand';                                           % DE for exploring agents should pull towards the element with the best objective function value
opt.social_DE_strategy ='DE/current-to-rand/1';                             % DE for social agents
opt.explore_all = 1;                                                        % all agents should perform local search
opt.v = 0;
opt.int_arch_mult=1;
opt.dyn_pat_search = 1;
opt.upd_subproblems = 0;
opt.max_rho_contr = 5;
opt.pat_search_strategy = 'standard';
opt.optimal_control = 1;
opt.vars_to_opt = [1;control_vars];
opt.oc.structure = structure;
opt.oc.smooth_scal_constr_fun = smooth_scal_const;
opt.oc.init_type = 'copy_ic';
opt.oc.x_0 = x_0;
opt.oc.x_f = x_f;
opt.oc.imposed_final_states = imposed_final_states;
opt.oc.state_vars = [0;state_vars];
opt.oc.control_vars = [0;control_vars];
%opt.oc.control_vars(1) = 0;


tol_conv = 1e-6;
maxits = 10000;
fminconoptions = optimset('Display','off','MaxFunEvals',maxits,'MaxIter',maxits,'TolCon',tol_conv,'TolX',1e-15,'GradConstr','off','Algorithm','sqp','MaxSQPIter', 10*length(vlb));


%% OPTIMISATION LOOP

for i=1:1
    
    mem(i).memory=macs7v16OC(@Simple_model_MACS_MOO,[],vlb,vub,opt,[],[],vlb,vub,structure,x_0,x_f,fminconoptions);    
    
end

%denormalise

%% plot 

[~,b] = sort(mem(1).memory(:,length(vlb)+1));
 
qq = mem(1).memory(b,:);    %sort wrt t_f

for i = 1:size(qq,1)
    
    xx = qq(i,1:length(vlb));
    %xx = xx.*structure.scale_optimisation_vars+structure.offset_optimisation_vars;
    [x,u,xb] = extract_solution(xx(2:end),structure,x_f);
    plot_solution_vs_time(x,u,x_0,xb,t_0,xx(1),structure.uniform_els,structure,i+1)
    %subplot(2,1,1)
    %axis([0 250 -pi/2 pi/2])
    %subplot(2,1,2)
    %axis([0 250 -1 2])
    drawnow
    
end
