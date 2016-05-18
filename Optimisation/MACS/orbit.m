function dx = orbit(t,x,mu)

%x(1) = x, x(2) = y, x(3) = z, x(4) = xdot, x(5) = ydot, x(6) = zdot

%mu = 3.9860e+14;

dx = zeros(6,1);

r = [x(1);x(2);x(3)];   %position vector
magr = norm(r);
%v = [x(4);x(5);x(6)];   % velocity vector
%magv = norm(v);

dx(1) = x(4);
dx(2) = x(5);
dx(3) = x(6);

dx(4) = -mu/magr^3*x(1);
dx(5) = -mu/magr^3*x(2);
dx(6) = -mu/magr^3*x(3);


end

