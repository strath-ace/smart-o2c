%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example of run of Multi-objective optimal control with MACS: Ascent
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

% Add path to optimiser folder
if isunix
    addpath(genpath('../../Optimisation'))
    addpath(genpath('../../Transcription'))
    addpath(genpath('MO_control/Ascent'))
else
    addpath(genpath('..\..\Optimisation'))
    addpath(genpath('..\..\Transcription'))    
    addpath(genpath('MO_control\Ascent'))
end

%% DFET PROBLEM TRANSCRIPTION

% Dynamics

f = @ascent_state_equation;
dfx = [];
dfu = [];
dft = [];

% Objective functions and Bolza's problem weights

g = @(x,u,t) [t 0; -x(2) 0];
weights = [1 0; 1 0];

% path constraints
c = [];

% initial and final conditions
t_0 = 0;

x_0 = [0 0 0 0];

imposed_final_states = [0 0 1 1];
x_f = [0 0.1 10 0];

%% Discretisation settings

num_elems = 4;
state_order = 6;
control_order = 6;
DFET = 1;
state_distrib = 'Lobatto';
control_distrib = 'Legendre';

integr_order = 2*state_order;

num_eqs = length(x_0);
num_controls = 1;

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

structure.f = f;
structure.dfx = dfx;
structure.dfu = dfu;
structure.dft = [];
structure.g = g;
structure.weights = weights;


state_bounds = [-inf inf;-inf inf; -inf inf;-inf inf];
control_bounds = [-pi/2 pi/2];

[vlb,vub] = transcribe_bounds(state_bounds,control_bounds,structure);

vlb = [100;vlb]';
vub = [250;vub]';

tol_conv = 1e-6;
maxits = 10000;

fminconoptions = optimoptions('fmincon','Display','off','MaxFunEvals',maxits,'TolCon',tol_conv,'GradConstr','on','Algorithm','interior-point','AlwaysHonorConstraints','none');

%% MACS PARAMETERS

opt.maxnfeval=10000;                                                        % maximum number of f evals 
opt.popsize=10;                                                             % popsize (for each archive)
opt.rhoini=1;                                                               % initial span of each local hypercube (1=full domain)
opt.F=0.9;                                                                  % F, the parameter for Differential Evolution
opt.CR=0.9;                                                                 % CR, crossover probability
opt.p_social=1;                                                             % popratio
opt.max_arch=10;                                                            % archive size
opt.coord_ratio=1;                                                          % ratio of total coordinates on which pattern search will be performed
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

% optimal control specific settings
opt.optimal_control = 1;
opt.vars_to_opt = ~isinf(vub);
opt.oc.structure = structure;
opt.oc.smooth_scal_constr_fun = [];
opt.oc.init_type = 'copy_ic';
opt.oc.x_0 = x_0;
opt.oc.x_f = x_f;
opt.oc.imposed_final_states = imposed_final_states;
opt.oc.state_vars = isinf(vub);
opt.oc.control_vars = ~isinf(vub);
opt.oc.control_vars(1) = 0;

%% OPTIMISATION LOOP

[x,fval,exitflag,output] = optimise_macs(@ascent_MOO,vlb,vub,opt,vlb,vub,structure,x_0,x_f,fminconoptions); 

%% plot 

[~,b] = sort(output.memory(:,length(vlb)+1));
 
qq = output.memory(b,:);    %sort wrt t_f

for i = 1:size(qq,1)
    
    [x,u,xb] = extract_solution(qq(i,2:length(vlb)),structure,x_f);
    plot_solution_vs_time(x,u,x_0,xb,t_0,qq(i,1),structure.uniform_els,structure,i+1)
    subplot(2,1,1)
    axis([0 250 0 120])
    subplot(2,1,2)
    axis([0 250 -pi pi])
    drawnow
    
end
