% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%

%% Phase 1

constants.Tmax = 193.044;
constants.g = 32.174;
constants.sigma = 5.49153484923381010e-5;
constants.c = 1580.9425279876559;
constants.h0 = 23800;

%% Definition of bounds on states, controls, times (i.e. transcription variables)

state_bounds = [-1 1e6; -1e4 1e4; 0.5 4];      % h, v, m
control_bounds = [];                         % 
time_bounds = [0 45];                        % this is the maximum range of time, t0 and tf must be within this range

%% Definition of initial conditions, final conditions, static variables, with their bounds if free

imposed_t_0 = 1;                                % t_0 is fixed
t_0 = 0;                                        % if t_0 is fixed, impose this value, otherwise it can used as initial guess
t0_bounds = [];                                 % if empty, use time_bounds. Not used if t_0 is imposed

imposed_t_f = 0;                                % t_f is free
t_f = 10;                                       % if t_f is fixed, impose this value, otherwise it can be used as initial guess
tf_bounds = [1 45];                             % if empty, use time_bounds. Not used if t_f is imposed

imposed_initial_states = ones(3,1);             % mask vector with the variables on which an initial condition is imposed
x_0 = [0 0 3]';                                 % as for t0, tf
x0_bounds = [];                                 % if empty, use state_bounds

imposed_final_states = zeros(3,1);              % mask vector with the variables on which a final condition is imposed
x_f = [0 0 0]';                                 % as for t0, tf
xf_bounds = [];                

other_vars_bounds = [];

%% Definition of initial guesses for static parameters and controls

other_vars_guess = [];
control_guess = [];

%% DFET Discretisation settings

num_elems = 3;
state_order = 9;
control_order = 0;
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

f = 'Goddard2_phase1_state_equations';
dfx = [];   % auto-compute :)
dfu = [];   % auto-compute :)

%% Definition of objective functions and Bolza's problem weights

g = [];
weights = [];%[0 0];
dgu0 = [];   % auto-compute :)
dgxf = [];   % auto-compute :)
dguf = [];   % auto-compute :)
dgxi = [];   % auto-compute :)
dgui = [];   % auto-compute :)

%% Definition of constraints

% Inequality path constraints

c = [];
dcx = [];   % auto-compute :)
dcu = [];   % auto-compute :)

% Equality path constraints

e = [];
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

q = [];
dqu0 = [];% auto-compute :)
dqxf = [];% auto-compute :)
dquf = [];% auto-compute :)
dqxi = [];% auto-compute :)
dqui = [];% auto-compute :)
wq = [];

%% Definition of type of initial guess (constant-linear interpolation/DFET)

init_type = 'DFET-all';

%% Definition of next phase(s) and the respective linkage functions

next_phases = [2];

link_funcs{1} = @Goddard2_link_ph_1_2;