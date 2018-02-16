% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%

%% Phase 4

constants.mu = 0.1407645794e17;
constants.Tc = 1.2;
constants.Isp = 300;
constants.g0 = 32.174;

%% Definition of bounds on states, controls, times (i.e. transcription variables)

p4 = 138334442.2575590;

pl = 3771261.588516088;
pu = 691672211.2877948;
Ll = 460*pi/180;
Lu = 641*pi/180;

state_bounds = [pl pu ; -1 1; -1 1; -1 1; -1 1; Ll Lu; 0.01 1.1];           % p, f, g, h, k, L
control_bounds = [-10*pi/180 10*pi/180; -90*pi/180 90*pi/180];                                            % theta, psi
time_bounds = [0 1e5];                                                      % this is the maximum range of time, t0 and tf must be within this range

%% Definition of initial conditions, final conditions, static variables, with their bounds if free

imposed_t_0 = 0;                                % t_0 is fixed
t_0 = 1e4;                                        % if t_0 is fixed, impose this value, otherwise it can used as initial guess
t0_bounds = [];                                 % if empty, use time_bounds. Not used if t_0 is imposed

imposed_t_f = 0;                                % t_f is free
t_f = 1.01e4;                                       % if t_f is fixed, impose this value, otherwise it can be used as initial guess
tf_bounds = [];                                 % if empty, use time_bounds. Not used if t_f is imposed

imposed_initial_states = [0 0 0 0 0 0 0];             % mask vector with the variables on which an initial condition is imposed
x_0 = [pl+1 0 0 0 0 Ll+1 1]';                       % as for t0, tf
x0_bounds = [];                                 % if empty, use state_bounds

imposed_final_states = [1 1 1 1 1 0 0];              % mask vector with the variables on which a final condition is imposed
x_f = [p4 0 0 0 0 0 0]';                               % as for t0, tf
xf_bounds = [];                

other_vars_bounds = [];

%% Definition of initial guesses for static parameters and controls

other_vars_guess = [];
control_guess = [0 0]';

%% DFET Discretisation settings

num_elems = 1;
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

f = 'Two_burns_phase2_state_equations';
dfx = [];   % auto-compute :)
dfu = [];   % auto-compute :)

%% Definition of objective functions and Bolza's problem weights

g = 'Two_burns_phase4_objective';
weights = [1 0];
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

h = 'Two_burns_bi_constraints';
dhu0 = [];% auto-compute :)
dhxf = [];% auto-compute :)
dhuf = [];% auto-compute :)
dhxi = [];% auto-compute :)
dhui = [];% auto-compute :)
wh = [1 0];

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