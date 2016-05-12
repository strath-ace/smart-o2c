function M = ascent_drag_mass_control_jacobian(x,u,t)

g = 0.0016; 
T = 4e-3;
rho0 = 0.1;
S = 1;
Cd = 1;
c = 1.5;

M = [0; -T*sin(u(1)); 0; T*cos(u(1)); 0];

end