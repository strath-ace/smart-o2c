% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%

%% Phase 1

% Formulation from Betts

%% Definition of problem specific constants

constants.m = 203000/32.174;
constants.rho0 =  0.002378;
constants.hr = 23800;
constants.Re = 20902900;
constants.S = 2690;
constants.a0 = -0.20704;
constants.a1 = 0.029244;
constants.mu = 0.14076539e17;
constants.b0 = 0.07854;
constants.b1 = -0.61592e-2;
constants.b2 = 0.621408e-3;
constants.c0 = 1.0672181;
constants.c1 = -0.19213774e-1;
constants.c2 = 0.21286289e-3;
constants.c3 = -0.10117249e-5;
constants.zeta = 0.2*pi/180;    % max rate of change for gamma and psi

% BC

hi = 260000;            % [ft]
phii = 0;               % [rad]
thetai = 0;             % [rad]
vi = 25600;             % [ft/s]
gammai = -1*pi/180;     % [rad]
psii = 90*pi/180;       % [rad]

hf = 80000;             % [ft]
thetaf = 15*pi/180;     % []
vf = 2500;              % [rad/s]
gammaf = -5*pi/180;     % [rad]

%% Definition of bounds on states, controls, times (i.e. transcription variables)

state_bounds = [0 400000; -179*pi/180 179*pi/180; -89*pi/180 89*pi/180; 1 30000; -89*pi/180 89*pi/180; -179*pi/180 179*pi/180];     % h [ft], lon (phi) [rad], lat (theta) [rad], v [ft/s], fpa (gamma) [rad], bank (psi) [rad]
control_bounds = [-90*pi/180 90*pi/180; -89*pi/180 1*pi/180];               % alpha [rad], beta [rad]
time_bounds = [0 2500];                        % this is the maximum range of time, t0 and tf must be within this range

%% Definition of initial conditions, final conditions, static variables, with their bounds if free

imposed_t_0 = 1;                                % t_0 is fixed
t_0 = 0;                                        % if t_0 is fixed, impose this value, otherwise it is ignored
t0_bounds = [];                                 % if empty, use time_bounds. Not used if t_0 is imposed

imposed_t_f = 0;                                % t_f is free
t_f = 2000;                                       % if t_f is fixed, impose this value, otherwise it is ignored
tf_bounds = [10 2500];                             % if empty, use time_bounds. Not used if t_f is imposed

imposed_initial_states = [1 1 1 1 1 1];%ones(6,1);           % mask vector with the variables on which an initial condition is imposed
x_0 = [hi phii thetai vi gammai psii]';                         
x0_bounds = [];                                 % if empty, use state_bounds

imposed_final_states = [1 0 0 1 1 0]';             % mask vector with the variables on which a final condition is imposed
x_f = [hf 0 0 vf gammaf 0]';                
xf_bounds = [0 400000; -179*pi/180 179*pi/180; -1*pi/180 89*pi/180; 1 30000; -89*pi/180 89*pi/180; -179*pi/180 179*pi/180];                

other_vars_bounds = [0 70];                   % max peak heat flux
other_vars_guess = [30];
control_guess = [0 0]';                                       

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

f = 'betts_descent_state_equations';
dfx = [];   % auto-compute :)
dfu = [];   % auto-compute :)

%% Definition of objective functions and Bolza's problem weights (as before)

g = 'betts_descent_objective';
weights = [1 0; 1 0];
dgu0 = [];   % auto-compute :)
dgxf = [];   % auto-compute :)
dguf = [];   % auto-compute :)
dgxi = [];   % auto-compute :)
dgui = [];   % auto-compute :)

%% Definition of constraints

% Inequality path constraints

c = 'betts_descent_ipath_constraints';
%c = @(x,u,t,static,scales) betts_descent_ipath_constraints(x,u,t,static,scales,constants);
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

