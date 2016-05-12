function M = state_equation(x,u,t)

g = 0.0016; 
T = 4e-3;
rho0 = 0.1;
S = 1;
Cd = 1;
c = 1.5;

M = [x(2); (T*cos((u(1)))-0.5*rho0*exp(-x(3))*S*Cd*(x(2).^2+x(4).^2).^(0.5)*x(2))/x(5); x(4); (T*sin((u(1)))-0.5*rho0*exp(-x(3))*S*Cd*(x(2).^2+x(4).^2).^(0.5)*x(4))/x(5)-g; -T/c];


end