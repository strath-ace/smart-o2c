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

constants.n = 0.06511*pi/180;    % orbital rotational rate [rad/s]
constants.hmax = 1e4;    % max angular momentum provided by cmg [ft*lbf*sec]
constants.umax = 200;
constants.J = [ 2.80701911616e7 4.822509936e5 -1.71675094448e7
                4.822509936e5 9.5144639344e7 6.02604448e4
                -1.71675094448e7 6.02604448e4 7.6594401336e7]; % ISS Inertia tensor [slug*ft^2]
constants.Jinv = inv(constants.J);  % avoid inverting it again, it's a constant!

constants.vcross = @(v) [0 -v(3) v(2);
                         v(3) 0 -v(1);
                         -v(2) v(1) 0];

% IC
r1i = 2.9963689649816e-3;
r2i = 1.5334477761054e-1;
r3i = 3.8359805613992e-3;

omega1i = -9.5380685844896e-6;
omega2i = -1.1363312657036e-3;
omega3i = 5.3472801108427e-6;

h1i = 5000;
h2i = 5000;
h3i = 5000;

% EC

h1f = 0;
h2f = 0;
h3f = 0;

%% Definition of bounds on states, controls, times (i.e. transcription variables)

state_bounds = [-1 1; -1 1; -1 1; -1e-2 1e-2; -1e-2 1e-2; -1e-2 1e-2; -1e4 1e4; -1e4 1e4; -1e4 1e4];     % q1, q2, q3, q4, w1, w2, w3, h1, h2, h3
control_bounds = [-200 200; -200 200; -200 200];               % ux, uy, uz
time_bounds = [0 1800];                        % this is the maximum range of time, t0 and tf must be within this range

%% Definition of initial conditions, final conditions, static variables, with their bounds if free

imposed_t_0 = 1;                                % t_0 is fixed
t_0 = 0;                                        % if t_0 is fixed, impose this value, otherwise it is ignored
t0_bounds = [];                                 % if empty, use time_bounds. Not used if t_0 is imposed

imposed_t_f = 1;                                % t_f is free
t_f = 1800;                                       % if t_f is fixed, impose this value, otherwise it is ignored
tf_bounds = [];                             % if empty, use time_bounds. Not used if t_f is imposed

imposed_initial_states = ones(9,1);           % mask vector with the variables on which an initial condition is imposed
x_0 = [r1i r2i r3i omega1i omega2i omega3i h1i h2i h3i]';                         
x0_bounds = [];                                 % if empty, use state_bounds

imposed_final_states = [0 0 0 0 0 0 1 1 1]';             % mask vector with the variables on which a final condition is imposed
x_f = [0 0 0 0 0 0 h1f h2f h3f]';                
xf_bounds = [];                

other_vars_bounds = [0 constants.hmax];       % bounds on gamma (max actuator effort)
other_vars_guess = [constants.hmax-1];
control_guess = [0 0 0]';

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

f = 'ZPM_alt_state_equations';
dfx = [];   % auto-compute :)
dfu = [];   % auto-compute :)

%% Definition of objective functions and Bolza's problem weights (as before)

g = 'ZPM_alt_objective';
weights = [1 0; 0 1];
dgu0 = [];   % auto-compute :)
dgxf = [];   % auto-compute :)
dguf = [];   % auto-compute :)
dgxi = [];   % auto-compute :)
dgui = [];   % auto-compute :)

%% Definition of constraints

% Inequality path constraints

c = 'ZPM_alt_ipath_constraints';
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

q = 'ZPM_alt_eepath_constraints';
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

