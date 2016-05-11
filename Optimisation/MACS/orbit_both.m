function dx = orbit_both(t,x,mt,mc,mu,erad,Ac,At,cdc,cdt,T,theta,phi)

%x(1) = xchaser, x(2) = ychaser, x(3) = zchaser
%x(4) = xtarget, x(5) = ytarget, x(6) = ztarget
%x(7) = xdotchaser, x(8) = ydotchaser, x(9) = zdotchaser
%x(10) = xdottarget, x(11) = ydottarget, x(12) = zdottarget

% theta = latitude, phi = longitude, as for initial conditions

%% Initialization

% if (norm(x(1:3))-erad)<0 || (norm(x(4:6))-erad)<0
%    keyboard 
% end

rho_c = exponential_atm_model((norm(x(1:3))-erad)*1e-3);
rho_t = exponential_atm_model((norm(x(4:6))-erad)*1e-3);

dx = zeros(12,1);

if length(x)<12
   keyboard 
end
rc = [x(1);x(2);x(3)];   %position vector of chaser
rt = [x(4);x(5);x(6)];   %position vector of target

magrc = norm(rc);
magrt = norm(rt);

vc = [x(7);x(8);x(9)];   % velocity vector of chaser
vt = [x(10);x(11);x(12)];   % velocity vector of chaser

magvc = norm(vc);
magvt = norm(vt);

%% Thrust calculation

xthrust = T*cos(theta)*cos(phi);
ythrust = T*cos(theta)*sin(phi);
zthrust = T*sin(theta);

%% Aerodynamic forces

aero_chaser_x = 0;
aero_chaser_y = 0;
aero_chaser_z = 0;

if magvc>0
   
    aero_chaser_x = -0.5*rho_c*Ac*cdc*magvc^2/magvc*x(7);
    aero_chaser_y = -0.5*rho_c*Ac*cdc*magvc^2/magvc*x(8);
    aero_chaser_z = -0.5*rho_c*Ac*cdc*magvc^2/magvc*x(9);
    
end

aero_target_x = 0;
aero_target_y = 0;
aero_target_z = 0;

if magvt>0
  
    aero_target_x = -0.5*rho_t*At*cdt*magvt^2/magvt*x(10);
    aero_target_y = -0.5*rho_t*At*cdt*magvt^2/magvt*x(11);
    aero_target_z = -0.5*rho_t*At*cdt*magvt^2/magvt*x(12);
    
end

%% Gravity forces

grav_chaser_x = -mu/magrc^3*x(1);
grav_chaser_y = -mu/magrc^3*x(2);
grav_chaser_z = -mu/magrc^3*x(3);

grav_target_x = -mu/magrt^3*x(4);
grav_target_y = -mu/magrt^3*x(5);
grav_target_z = -mu/magrt^3*x(6);

%% Computation of time derivatives

dx(1) = x(7);
dx(2) = x(8);
dx(3) = x(9);

dx(4) = x(10);
dx(5) = x(11);
dx(6) = x(12);

dx(7) = grav_chaser_x+(aero_chaser_x+xthrust)/mc;
dx(8) = grav_chaser_y+(aero_chaser_y+ythrust)/mc;
dx(9) = grav_chaser_z+(aero_chaser_z+zthrust)/mc;

dx(10) = grav_target_x+aero_target_x/mt;
dx(11) = grav_target_y+aero_target_y/mt;
dx(12) = grav_target_z+aero_target_z/mt;

end

