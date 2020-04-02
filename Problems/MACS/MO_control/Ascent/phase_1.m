% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%

%% Phase 1

% formulation from Betts

%% Definition of problem specific constants

constants.g = 1.6e-3;
constants.a = 4e-3;

% EC
yf = 10;
vyf = 0;

%% Definition of bounds on states, controls, times (i.e. transcription variables)

state_bounds = [-1 200; -10 10; -1 11; -10 10];     % r, vr, theta, vt
control_bounds = [-pi/2 pi/2];                       % u (rad)
time_bounds = [0 250];                        % this is the maximum range of time, t0 and tf must be within this range

%% Definition of initial conditions, final conditions, static variables, with their bounds if free

imposed_t_0 = 1;                                % t_0 is fixed
t_0 = 0;                                        % if t_0 is fixed, impose this value, otherwise it is ignored
t0_bounds = [];                                 % if empty, use time_bounds. Not used if t_0 is imposed

imposed_t_f = 0;                                % t_f is free
t_f = 150;                                       % if t_f is fixed, impose this value, otherwise it is ignored
tf_bounds = [100 250];                            % if empty, use time_bounds. Not used if t_f is imposed

imposed_initial_states = [1 1 1 1];           % mask vector with the variables on which an initial condition is imposed
x_0 = [0 0 0 0]';                         
x0_bounds = [];                                 % if empty, use state_bounds

imposed_final_states = [0 0 1 1];             % mask vector with the variables on which a final condition is imposed
x_f = [0 0 yf vyf]';                
xf_bounds = [];                

other_vars_bounds = [];                         % no static vars
other_vars_guess = [];
control_guess = 0;                                       

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

f = 'ascent_state_equations';
%f = @(x,u,t,static,scales) ascent_state_equations(x,u,t,static,scales,constants);
dfx = [];   % auto-compute :)
dfu = [];   % auto-compute :)

%% Definition of objective functions and Bolza's problem weights (as before)

g = 'ascent_objective';
weights = [1 0; 1 0];
dgu0 = [];   % auto-compute :)
dgxf = [];   % auto-compute :)
dguf = [];   % auto-compute :)
dgxi = [];   % auto-compute :)
dgui = [];   % auto-compute :)

%% Definition of constraints

% Inequality path constraints

c = [];%
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

next_phases = [];
