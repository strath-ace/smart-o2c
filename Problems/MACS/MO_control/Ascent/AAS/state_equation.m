function dx = state_equation(x,u,t)

% initialisation
% states:   x(1) = h, altitude
%           x(2) = theta, longitude
%           x(3) = v, modulus of velocity
%           x(4) = gamma, flight path angle
%           x(5) = m, total mass
%
% controls: u(1) = delta, throttle
%           u(2) = alpha, angle of attack in degrees (BEWARE: convert to radians before using them in sin and cos functions!!!)

dx = zeros(5,1);

% de-normalise input



% constants

r_e = 6371e3;                       % Earth radius, [m]
mdot_fuel = 850;                    % Max fuel flow, [kg/s]
g_0 = 9.81;                         % Gravitational acceleration on Earth's surface [m/s^2]
I_sp = 450;                         % Specific impulse [s]
rho_20km = 8.80349E-2;              % Air density at 20km, from U.S. Standard Atmosphere, 1976 [kg/m^3]
rho_40km = 3.85101E-3;              % Air density at 40km, from U.S. Standard Atmosphere, 1976 [kg/m^3]
c = -log(rho_40km/rho_20km)/2e4;    % Coefficient for exponential fit between the two density values [1/m]
C_d_0 = 0.012;                      % Zero lift drag coefficient []
k = 2.09;                           % Oswald factor []
C_l_alpha = 0.012;                  % Slope of CL(alpha) curve []
omega_e = 7.2921150e-5;             % Earth angular velocity wrt its axis [rad/s]
S = 200;                            % Surface area [m^2]
mu_e = 398600e9;                    % Earth gravitational constant [m^3/s^2]

% computation of forces, densities, etc... tidier to perform it separately from the dynamics

rho = rho_20km*exp(-(x(1)-2e4)*c);                  % Air density, BEWARE: MUST GIVE rho_20km at h=20km, thus the difference (x(1)-2e4)
F_T = mdot_fuel*g_0*I_sp*u(1);                      % Thrust
L = 0.5*rho*x(3)^2*S*C_l_alpha*u(2);                % Lift
D = 0.5*rho*x(3)^2*S*(C_d_0+k*(C_l_alpha*u(2))^2);  % Drag

% the actual dynamics

dx(1) = x(3)*sin(x(4));
dx(2) = x(3)*cos(x(4))/(r_e+x(1));
dx(3) = (F_T*cos(pi/180*u(2))-D)/x(5) - (mu_e/(x(1)+r_e)^2 - omega_e^2*(x(1)+r_e))*sin(x(4));
dx(4) = (F_T*sin(pi/180*u(2))+L)/(x(5)*x(3)) - (mu_e/((x(1)+r_e)^2*x(3))-x(3)/(x(1)+r_e) - omega_e^2*(x(1)+r_e)/x(3))*cos(x(4)) + 2*omega_e;
dx(5) = -mdot_fuel*u(1);

end