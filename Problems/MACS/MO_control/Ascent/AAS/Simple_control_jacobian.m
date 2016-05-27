function M = Simple_control_jacobian(x,u,t)

% initialisation
% states:   x(1) = h, altitude
%           x(2) = theta, longitude
%           x(3) = v, modulus of velocity
%           x(4) = gamma, flight path angle
%           x(5) = m, total mass
%
% controls: u(1) = delta, throttle
%           u(2) = alpha, angle of attack in degrees (BEWARE: convert to radians before using them in sin and cos functions!!!)

M = zeros(5,2);

% de normalising manually input :(

tscale = 570;
toffs = 0;

xscale = [50000 2*pi 999 pi 250000]';
xoffs = [0 0 1 -pi/2 50000]';

uscale = [1 60]';
uoffs = [0 -10]';

x = x.*xscale+xoffs;
u = u.*uscale+uoffs;

t = t*tscale+toffs;

% constants

%r_e = 6371e3;                       % Earth radius, [m]
mdot_fuel = 850;                    % Max fuel flow, [kg/s]
g_0 = 9.81;                         % Gravitational acceleration on Earth's surface [m/s^2]
I_sp = 450;                         % Specific impulse [s]
rho_20km = 8.80349E-2;              % Air density at 20km, from U.S. Standard Atmosphere, 1976 [kg/m^3]
rho_40km = 3.85101E-3;              % Air density at 40km, from U.S. Standard Atmosphere, 1976 [kg/m^3]
c = -log(rho_40km/rho_20km)/2e4;    % Coefficient for exponential fit between the two density values [1/m]
%C_d_0 = 0.012;                      % Zero lift drag coefficient []
k = 2.09;                           % Oswald factor []
C_l_alpha = 0.012;                  % Slope of CL(alpha) curve []
%omega_e = 7.2921150e-5;             % Earth angular velocity wrt its axis [rad/s]
S = 200;                            % Surface area [m^2]
%mu_e = 398600e9;                    % Earth gravitational constant [m^3/s^2]

% computation of forces, densities, etc... tidier to perform it separately
% from the Jacobian

rho = rho_20km*exp(-(x(1)-2e4)*c);                  % Air density, BEWARE: MUST GIVE rho_20km at h=20km, thus the difference (x(1)-2e4)
F_T = mdot_fuel*g_0*I_sp*u(1);                      % Thrust
%L = 0.5*rho*x(3)^2*S*C_l_alpha*u(2);                % Lift
%D = 0.5*rho*x(3)^2*S*(C_d_0+k*(C_l_alpha*u(2))^2);  % Drag

% computation of derivatives of forces, densities, etc... tidier to perform
% it separately from the Jacobian

dF_T_ddelta = mdot_fuel*g_0*I_sp;            % derivative of Thrust wrt throttle
dL_dalpha = 0.5*rho*x(3)^2*S*C_l_alpha;      % derivative of Lift wrt alpha
dD_dalpha = rho*x(3)^2*S*k*C_l_alpha^2*u(2); % derivative of Drag wrt alpha

% the actual Jacobian (has few nonzero elements)

M(3,1) = dF_T_ddelta*cos(pi/180*u(2))/x(5);                     % derivative of third eqn (velocity) wrt throttle
M(3,2) = (-F_T*pi/180*sin(pi/180*u(2)) - dD_dalpha)/x(5);       % derivative of third eqn (velocity) wrt angle of attack
M(4,1) = dF_T_ddelta*sin(pi/180*u(2))/(x(5)*x(3));              % derivative of fourth eqn (gamma) wrt throttle
M(4,2) = (F_T*pi/180*cos(pi/180*u(2)) + dL_dalpha)/(x(5)*x(3)); % derivative of fourth eqn (gamma) wrt angle of attack
M(5,1) = -mdot_fuel;                                            % derivative of fifth eqn (mass) wrt throttle

% renormalise output

scalemat = repmat(uscale',size(xscale,1),1)./repmat(xscale,1,size(uscale,1));

M = M.*scalemat;
end