% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%

%% Phase 1

%% Definition of problem specific constants

constants.a = 7071000;
constants.mu = 398e12;
constants.n = sqrt(constants.mu/constants.a^3);
constants.m = 100;
constants.alphamax = 0.1;
constants.taumax = 1;
constants.J11 = 1000;
constants.J22 = 2000;
constants.J33 = 1000;
constants.K11 = 1000;
constants.K22 = 2000;
constants.K33 = 1000;
constants.avec = [0 1.01 0]';
constants.bvec = [0 -1.01 0]';

%% Definition of bounds on states, controls, times (i.e. transcription variables)

state_bounds = [-20 20; -20 20; -20 20; -20 20; -20 20; -20 20; -1.1 1.1; -1.1 1.1; -1.1 1.1; -1.1 1.1; -20*pi 20*pi; -20*pi 20*pi; -20*pi 20*pi; -1.1 1.1; -1.1 1.1; -1.1 1.1; -1.1 1.1; -20*pi 20*pi; -20*pi 20*pi; -20*pi 20*pi]; % x, y, z, vx, vy, vz, q1, q2, q3, q4, w1, w2, w3, p1, p2, p3, p4, phi1, phi2, phi3
control_bounds = [-1 1; -1 1; -1 1; -1 1; -1 1; -1 1];          % alpha, tau
time_bounds = [0 420];                        % this is the maximum range of time, t0 and tf must be within this range

%% Definition of initial conditions, final conditions, static variables, with their bounds if free

imposed_t_0 = 1;                                % t_0 is fixed
t_0 = 0;                                        % if t_0 is fixed, impose this value, otherwise it is ignored
t0_bounds = [];                                 % if empty, use time_bounds. Not used if t_0 is imposed

imposed_t_f = 0;                                % t_f is free
t_f = 100;                                       % if t_f is fixed, impose this value, otherwise it is ignored
tf_bounds = [1 420];                            % if empty, use time_bounds. Not used if t_f is imposed

imposed_initial_states = ones(20,1);           % mask vector with the variables on which an initial condition is imposed
x_0 = [0 -10 0 0 0 0 0 0 0 1 0 0 0 -0.05 0 0 (1-0.05^2)^0.5 0 0.0349 0.017453]';                         
x0_bounds = [];                                 % if empty, use state_bounds

imposed_final_states = zeros(20,1);             % mask vector with the variables on which a final condition is imposed
x_f = x_0;                
xf_bounds = [];                

other_vars_bounds = [];                         % no static vars
other_vars_guess = [];
control_guess = [-0.95 0 0.95 0 0 0]';                                       

%% DFET Discretisation settings

num_elems = 4;
state_order = 6;
control_order = 6;
DFET = 1;
state_distrib = 'Bernstein';
control_distrib = 'Bernstein';
test_distrib = 'Bernstein';
integr_type = 'Legendre';
num_eqs = size(state_bounds,1);
num_controls = size(control_bounds,1);

%% Definition of dynamics (can be moved after the initial transcription, so the scaling factors are computed internally and the user cannot mess them up, and they are ensured to be consistent through all the code)

% Wrapping commonly used constants into a handy and clean structure. 
% Could also include functions/models, like atmospheric models, etc
% Very useful and flexible.

f = 'Docking_state_equations';
dfx = [];   % auto-compute :)
dfu = [];   % auto-compute :)

%% Definition of objective functions and Bolza's problem weights (as before)

g = 'Docking_objective';
weights = [0 1];
dgu0 = [];   % auto-compute :)
dgxf = [];   % auto-compute :)
dguf = [];   % auto-compute :)
dgxi = [];   % auto-compute :)
dgui = [];   % auto-compute :)


%% Definition of constraints

% Inequality path constraints

c = 'Docking_ipath_constraints';
dcx = [];   % auto-compute :)
dcu = [];   % auto-compute :)

% Equality path constraints

e =   [];%@(x,u,t,static,scales) Docking_epath_constraints(x,u,t,static,scales,constants);
dex = [];   % auto-compute :)
deu = [];   % auto-compute :)

% Final time inequality constraints and integral inequality constraints

h = [];
dhu0 = [];% auto-compute :)
dhxf = [];% auto-compute :)
dhuf = [];% auto-compute :)
dhxi = [];% auto-compute :)
dhui = [];% auto-compute :)
wh = [];

% Final time equality constraints and integral equality contraints

q = 'Docking_bc';
dqu0 = [];% auto-compute :)
dqxf = [];% auto-compute :)
dquf = [];% auto-compute :)
dqxi = [];% auto-compute :)
dqui = [];% auto-compute :)
wq = [ones(6,1) zeros(6,1)];

%% Definition of type of initial guess (constant-linear interpolation/DFET)

init_type = 'DFET-all';

%% Definition of next phase(s) and the respective linkage functions

next_phases = [];
