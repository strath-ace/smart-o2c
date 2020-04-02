% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%

%% Phase 2

%% Definition of problem specific constants

%constants.m = 92e3;                    % Vehicle mass [kg]
constants.Re = 6371e3;                  % Earth radius [m]
constants.mu = 3.986004418e14;          % Earth's gravitational parameter [m^3/s^2]
%constants.S = 50;                      % Vehicle wing area [m^2]
constants.c0 = 1.0672181;               % Thermal coefficients, UNKNOWN UNITS, COMPUTATIONS WILL BE PERFORMED IN IMPERIAL SYSTEM AND CONVERTED BACK
constants.c1 = -0.19213774e-1;          % Thermal coefficients, UNKNOWN UNITS, COMPUTATIONS WILL BE PERFORMED IN IMPERIAL SYSTEM AND CONVERTED BACK
constants.c2 = 0.21286289e-3;           % Thermal coefficients, UNKNOWN UNITS, COMPUTATIONS WILL BE PERFORMED IN IMPERIAL SYSTEM AND CONVERTED BACK
constants.c3 = -0.10117249e-5;          % Thermal coefficients, UNKNOWN UNITS, COMPUTATIONS WILL BE PERFORMED IN IMPERIAL SYSTEM AND CONVERTED BACK
constants.qmax = 800e3;                 % Heat flux density threshold [w/m^2]
constants.zeta = 0.2*pi/180;            % max rate of change for gamma and psi [rad/s]
constants.maxacc = 3*9.81;              % max acceleration [m/s^2]

constants.maxThrust = 1340e3;           % Thrust of the engine at sea level
constants.Isp = 450;%336;                    % specific impulse of engine
constants.Ae = 1.9607;                  % Nozzle exit area
constants.g0 = constants.mu/constants.Re^2; % gravitational acceleration at sea level (computed for consistency)
constants.omega_e = 7.2921159e-5;       % Earth angular velocity [rad/s]

constants.a = constants.Re+100e3;       % target orbit has 100 km semimajor axis
constants.e = 0;                        % zero eccentricity
constants.i = 0;                        % and zero inclination

% autobuild function :D

load('X34aero/cl_model_values');
load('X34aero/cl_model_names');

for j = 1:length(cl_model_values)

    eval([cl_model_names{j},'=',num2str(cl_model_values(j)),';']);
    
end

cl_fun = @(x,y) (a0+a1*x+a2*x^2)+(b0+b1*x+b2*x^2)*(k1/l1)*((y-s1)/l1)^(k1-1)*exp(-((y-s1)/l1)^k1)+(c0+c1*x+c2*x^2)*(k2/l2)*((y-s3)/l2)^(k2-1)*exp(-((y-s3)/l2)^k2);

% autobuild function :D

load('X34aero/cd_model_values');
load('X34aero/cd_model_names');

for j = 1:length(cd_model_values)

    eval([cd_model_names{j},'=',num2str(cd_model_values(j)),';']);
    
end

cd_fun = @(x,y) (a0+a1*x+a2*x^2+a3*x^3)+(b0+b1*x+b2*x^2+b3*x^3)*(k1/l1)*((y-s1)/l1)^(k1-1)*exp(-((y-s1)/l1)^k1)+(c0+c1*x+c2*x^2+c3*x^3)*(k2/l2)*((y-s3)/l2)^(k2-1)*exp(-((y-s3)/l2)^k2);

constants.cl_fun = cl_fun;
constants.cd_fun = cd_fun;

% Initial conditions

hi = 10e3;             % [m]
phii = 0;               % [rad]
thetai = 0;             % [rad]
vi = 200;              % [m/s]
gammai = 0*pi/180;     % [rad]
psii = 90*pi/180;       % [rad]

%% Definition of bounds on states, controls, times (i.e. transcription variables)

state_bounds = [1 1.2e5; -1*pi/180 179*pi/180; -89*pi/180 89*pi/180; 1 1e4; -89*pi/180 89*pi/180; 0 pi; 10e3 200e3];     % h [m], lon (phi) [rad], lat (theta) [rad], v [m], fpa (gamma) [rad], bank (psi) [rad]
control_bounds = [0*pi/180 35*pi/180; -40*pi/180 40*pi/180; 0 1];           % alpha [rad], beta [rad], thrust
time_bounds = [0 2000];                                                         % this is the maximum range of time, t0 and tf must be within this range

%% Definition of initial conditions, final conditions, static variables, with their bounds if free

imposed_t_0 = 1;                                % t_0 is fixed
t_0 = 0;                                      % if t_0 is fixed, impose this value, otherwise it is ignored
t0_bounds = [];                          % if empty, use time_bounds. Not used if t_0 is imposed

imposed_t_f = 0;                                % t_f is free
t_f = 200;                                      % if t_f is fixed, impose this value, otherwise it is ignored
tf_bounds = [];                          % if empty, use time_bounds. Not used if t_f is imposed

imposed_initial_states = [0 0 0 0 0 0 0]';%ones(6,1);           % mask vector with the variables on which an initial condition is imposed
x_0 = [20e3 0 0 4e3 0 pi/2 40e3]';                         
x0_bounds = [];                                 % if empty, use state_bounds

imposed_final_states = [1 0 1 1 1 1 0]';             % mask vector with the variables on which a final condition is imposed
x_f = [100e3 0 0 sqrt(constants.mu/constants.a)-constants.omega_e*constants.Re 0 pi/2 0]';                   % 7.4km/s because there's also Earth'rotation (we're not in an Inertial Frame)
xf_bounds = [];

other_vars_bounds = [20 400; 10e3 200e3];                    % wing surface, phase 1 initial mass
other_vars_guess = [300;  70e3];
control_guess = [5*pi/180 0 0.9]';                                       

%% DFET Discretisation settings

num_elems = 3;
state_order = 7;
control_order = 7;
DFET = 1;
state_distrib = 'Bernstein';
control_distrib = 'Bernstein';
num_eqs = size(state_bounds,1);
num_controls = size(control_bounds,1);

%% Definition of dynamics (can be moved after the initial transcription, so the scaling factors are computed internally and the user cannot mess them up, and they are ensured to be consistent through all the code)

% Wrapping commonly used constants into a handy and clean structure. 
% Could also include functions/models, like atmospheric models, etc
% Very useful and flexible.

f = 'New4ESA_state_equations';
dfx = [];   % auto-compute :)
dfu = [];   % auto-compute :)

%% Definition of objective functions and Bolza's problem weights (as before)

g = [];%'New4ESA_objective';
weights = [];%[1/20 0];
dgu0 = [];   % auto-compute :)
dgxf = [];   % auto-compute :)
dguf = [];   % auto-compute :)
dgxi = [];   % auto-compute :)
dgui = [];   % auto-compute :)

%% Definition of constraints

% Inequality path constraints

c = 'New4ESA_ipath_constraints';
dcx = [];   % auto-compute :)
dcu = [];   % auto-compute :)

% Equality path constraints

e = [];
dex = [];   % auto-compute :)
deu = [];   % auto-compute :)

% Final time inequality constraints and integral inequality constraints

h = [];%'New4ESA_ei_constraints';
dhu0 = [];
dhxf = [];
dhuf = [];
dhxi = [];
dhui = [];
wh = [];

% Final time equality constraints and integral equality contraints

q = 'New4ESA_ee_constraints2';
dqu0 = [];
dqxf = [];
dquf = [];
dqxi = [];
dqui = [];
wq = [1 0];%[1 0; 1 0; 1 0];

%% Definition of type of initial guess (constant-linear interpolation/DFET)

init_type = 'DFET-all';

%% Definition of next phase(s) and the respective linkage functions

next_phases = [];
