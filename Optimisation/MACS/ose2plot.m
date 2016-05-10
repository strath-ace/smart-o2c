function ose2plot(in)

%% Simulation parameters
tmax = 7200;                                                                % [s]
toll = 1e3;                                                                 % [m]
ttoll = 1e-1;

%% Input decoding

xx = zeros(1,2);

% ttotal = in(1);
% delay = in(2);

delay = in(1)*tmax;
delay = delay*(delay>=ttoll)+0*(delay<ttoll);

vals = in(2:end);
n = length(vals)/4;
DT_on = vals(1:n)*tmax;
DT_on = DT_on.*(DT_on>=ttoll)+0*(DT_on<ttoll);
DT_off = vals(n+1:2*n)*tmax;
DT_off = DT_off.*(DT_off>=ttoll)+0*(DT_off<ttoll);
theta = vals(2*n+1:3*n)*2*pi;
phi = vals(3*n+1:4*n)*(pi/2);

Ttot = sum([delay DT_on DT_off]);

if Ttot>tmax
   
    xx(1) = realmax;
    xx(2) = realmax;
    c = realmax;
    return
end

T_on = zeros(1,n);
T_off= zeros(1,n);

T_on(1) = delay + DT_on(1);
T_off(1) = T_on(1)+DT_off(1);

if n>1
    for i = 2:n
        
        T_on(i) = T_off(i-1)+DT_on(i);
        T_off(i) = T_on(i)+DT_off(i);
        
    end
    
end

if T_off(end)>tmax
   keyboard 
end

%% Environmental parameters

rho = 1e-14;
grav_const = 398600*10^9;                                                   %[m3/s2]
earth_rad = 6371e3;                                                         %[m]
T = 10;                                                                      %[m/s^2], acceleration given by the thrusters

%% Initial condition of chaser and its physical parameters

% All positions and velocity are ABSOLUTE (wrt fixed earth)

chaser.altitude = 90e3;                                                     %[m, from surface]
chaser.latitude = 0*pi/180;                                                 %[rad, from equator]
chaser.longitude = 0*pi/180;                                                %[rad, from Greenwich]

r = chaser.altitude+earth_rad;
x = r*cos(chaser.latitude)*cos(chaser.longitude);
y = r*cos(chaser.latitude)*sin(chaser.longitude);
z = r*sin(chaser.latitude);

chaser.x0 = [x y z];                                                        %[m]
chaser.xdot0 = [0 (grav_const/r^3)^0.5*r 0];                                %[m/s]
chaser.mass = 1435;                                                         %[kg]
chaser.inertia = [2040 130 25; 130 1670 -55; 25 -55 2570];                  %[kg/m2]
chaser.bbox = [3;5;2];                                                      %[m]
chaser.cd = [2.9;2.9;2.9];                                                  %[]
chaser.cm = [0;0;0];                                                        %[];
chaser.area = [56.64;56.64;56.64];                                           %[m^2];

%% Initial condition of targer and its physical parameters

target.altitude = 200e3;                                                    %[m, from surface]
target.latitude = 0*pi/180;                                                 %[rad, from equator]
target.longitude = 0*pi/180;                                                %[rad, from Greenwich]

r = target.altitude+earth_rad;
x = r*cos(target.latitude)*cos(target.longitude);
y = r*cos(target.latitude)*sin(target.longitude);
z = r*sin(target.latitude);

target.x0 = [x y z];                                                        %[m]
target.xdot0 = [0 (grav_const/r^3)^0.5*r 0];                                %[m/s]
target.mass = 7792;                                                         %[kg]
target.inertia = [17012 401 -2167; 401 124725 345; -2167 345 129009];       %[kg/m2]
target.bbox = [4;5;27];                                                     %[m]
target.cd = [2.9;2.9;2.9];                                                  %[]
target.cm = [0;0;0];                                                        %[];
target.area = [56.64;56.64;56.64];                                          %[m^2];

% Simple orbiting of one body
% options = odeset ('Event',@(t,y) on_ground(t,y,earth_rad));
% [tchas,ychas] = ode45(@(t,y) orbit(t,y,grav_const),[0 7200],[chaser.x0;chaser.xdot0],options);
% [ttar,ytar] = ode45(@(t,y) orbit(t,y,grav_const),[0 7200],[target.x0;target.xdot0],options);

%% ODE solution on the arcs and assignment of objective values
clear x y z
options = odeset ('Event',@(t,y) on_ground_both(t,y,earth_rad,toll),'InitialStep',1e-1);

tout = 0;
yout(1,:) = [chaser.x0 target.x0 chaser.xdot0 target.xdot0];
tspan = [0 delay];
y0 = yout(1,:);

%initial arc with delay, if not null: engines off

if delay~=0
    
    % solve ode with engines off
    [t,y] = ode45(@(t,y) orbit_both(t,y,grav_const,rho,chaser.area(1),target.area(1),chaser.cd(1),target.cd(1),0,0,0),tspan,y0,options);
    
    % append solution
    tout = [tout; t(2:end)];
    yout = [yout; y(2:end,:)];
    
    % update temporary objective values
    xx(1) = tout(end);
    xx(2) = xx(2); %no fuel consumption in this case
    
    c = norm(yout(end,1:3)-yout(end,4:6))-toll;
    
    % failed mission
    if (norm(yout(end,1:3))<=earth_rad) || (norm(yout(end,4:6))<=earth_rad)
        xx(1) = realmax;
        xx(2) = realmax;
        c = realmax;
        return
    end
    
    %successful
    if norm(yout(end,1:3)-yout(end,4:6))<toll
        return
    end
    
end

tspan = [delay T_on(1)];
y0 = yout(end,:);

% n arcs, engines on and then off

for i=1:n
    
    if tspan(2)>tspan(1)
        % solve ode with engines on
        [t,y] = ode45(@(t,y) orbit_both(t,y,grav_const,rho,chaser.area(1),target.area(1),chaser.cd(1),target.cd(1),T,theta(i),phi(i)),tspan,y0,options);
        
        % append solution
        tout = [tout; t(2:end)];
        yout = [yout; y(2:end,:)];
        
        % compute temporary objective values
        xx(1) = tout(end);
        xx(2) = xx(2)+tout(end)-tspan(1);
        
        c = norm(yout(end,1:3)-yout(end,4:6))-toll;
        
        % failed mission
        if (norm(yout(end,1:3))<=earth_rad) || (norm(yout(end,4:6))<=earth_rad)
            xx(1) = realmax;
            xx(2) = realmax;            
            c = realmax;
            break
        end
        
        %successful
        if norm(yout(end,1:3)-yout(end,4:6))<toll
            break
        end
        
    end
    
    tspan = [T_on(i) T_off(i)];  %new time span
    y0 = yout(end,:);          %new initial condition
    
    if tspan(2)>tspan(1)
        
        % solve ode with engines off
        [t,y] = ode45(@(t,y) orbit_both(t,y,grav_const,rho,chaser.area(1),target.area(1),chaser.cd(1),target.cd(1),0,0,0),tspan,y0,options);
        
        % append solution
        tout = [tout; t(2:end)];
        yout = [yout; y(2:end,:)];
        
        % update temporary objective values
        xx(1) = tout(end);
        xx(2) = xx(2); %no fuel consumption in this case
        
        c = norm(yout(end,1:3)-yout(end,4:6))-toll;
        
        % failed mission
        if (norm(yout(end,1:3))<=earth_rad) || (norm(yout(end,4:6))<=earth_rad)
            xx(1) = realmax;
            xx(2) = realmax;            
            c = realmax;
            break
        end
        
        %successful
        if norm(yout(end,1:3)-yout(end,4:6))<toll
            break
        end
        
    end
    
    if i<n
        
        tspan = [T_off(i) T_on(i+1)];  %new time span
        y0 = yout(end,:);          %new initial condition
        
    end
    
end

% final arc with engines off, if times remains and we're not at position

% if tout(end)<tmax && (norm(yout(end,1:3))>earth_rad) && (norm(yout(end,4:6))>earth_rad) && ((norm(yout(end,1:3)-yout(end,4:6)))>toll)
%    
%     tspan = [tout(end) tmax];
%     y0 = yout(end,:);
%     
%     % solve ode with engines off
%     [t,y] = ode45(@(t,y) orbit_both(t,y,grav_const,rho,chaser.area(1),target.area(1),chaser.cd(1),target.cd(1),0,0,0),tspan,y0,options);
%     
%     % append solution
%     tout = [tout; t(2:end)];
%     yout = [yout; y(2:end,:)];
%     
%     % update temporary objective values
%     xx(1) = tout(end);
%     xx(2) = xx(2); %no fuel consumption in this case
%     
%     c = norm(yout(end,1:3)-yout(end,4:6))-toll;
%     
%     % failed mission
%     if (norm(yout(end,1:3))<=earth_rad) || (norm(yout(end,4:6))<=earth_rad)
%         xx(1) = realmax;
%         xx(2) = realmax;
%         c = realmax;
%         return
%     end
%     
%     %successful
%     if norm(yout(end,1:3)-yout(end,4:6))<toll
%         return
%     end
%     
% end


% degenerate solution, for some reason
if all(xx)==0
    
    xx(1) = realmax;
    xx(2) = realmax;
    %xx(3) = realmax;
    
    c = realmax;
    
end


% Nice 3D plot

[xsphere,ysphere,zsphere]=sphere(100);
xsphere = earth_rad*xsphere;
ysphere = earth_rad*ysphere;
zsphere = earth_rad*zsphere;
% 
% plot3(yout(:,1),yout(:,2),yout(:,3),'r','LineWidth',3);
% hold on
% plot3(yout(:,4),yout(:,5),yout(:,6),'g','LineWidth',3);
% axis equal
% surf(xsphere,ysphere,zsphere,ones(size(xsphere)))
% hold on
% plot3(yout(:,1),yout(:,2),yout(:,3),'g','LineWidth',3);
% plot3(yout(:,4),yout(:,5),yout(:,6),'r','LineWidth',3);
for i=1:size(tout,1)
    hold off
    plot3(yout(1:i,1),yout(1:i,2),yout(1:i,3),'g')%,'LineWidth',3);
    hold on
    surf(xsphere,ysphere,zsphere,ones(size(xsphere)))
    plot3(yout(1,1),yout(1,2),yout(1,3),'go')
    plot3(yout(end,1),yout(end,2),yout(end,3),'gx')
    axis equal
    plot3(yout(1:i,4),yout(1:i,5),yout(1:i,6),'r')%,'LineWidth',3);
    plot3(yout(1,4),yout(1,5),yout(1,6),'ro')
    plot3(yout(end,4),yout(end,5),yout(end,6),'rx')
    drawnow
    %surf(xsphere,ysphere,zsphere,ones(size(xsphere)))
    %legend('Chaser','','','Target','','')
end

%figure(2)
%plot(tout,(sum((yout(:,1:3)-yout(:,4:6)).^2,2)).^0.5);

end
