function theta = M2theta(M,e)
%Function M2theta
% Find true anomaly from mean anomaly and eccentricity
%
% INPUT
% M: mean anomaly
% e: eccentricity
%
% OUTPUT
% theta: true anomaly
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Niccolo' Gastaldello, November 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Computation

% Initial values
err = 1;                            %Set an initial, big error  
toll = 1.e-7;                       %Set the tolerance

if M < pi                           %Criterion for choosing initial guess 
    E = M + e/2;                    % (from Orb. Mech. for Eng. Students)
else
    E = M - e/2;
end

% Cycle
while err > toll;                  %Cycle until the difference is smaller
                                    %than tolerance
    ratio = (M - E + e*sin(E) ) / (1 - e*cos(E));
    E = E + ratio;
    err = abs(ratio);
end

% Obtain theta
tan_theta2 = sqrt((1+e) / (1-e))*tan(E/2);          %Calculate theta
theta = 2*atan(tan_theta2);
 
if theta < 0                         %Transform negative angle in positive
    theta = theta + 2*pi;
end
end

