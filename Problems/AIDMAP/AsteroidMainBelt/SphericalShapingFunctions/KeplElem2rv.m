function [r,v] = KeplElem2rv(a,e,i,w,Omega,theta,mu)

%KeplElem2rv: Function to obtain a cartesian state vector in the inertial frame
% starting from the Kepler orbital parameters
%
% INPUT
% Kepler parameters, mu (gravitational parameter)

% OUTPUT
% r,v: position and velocity

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Niccolo' Gastaldello, February 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Write the time derivatives
%------------------------------------------------------------------------%

p = a*(1-e^2);

R_pf = (p/(1+e*cos(theta)))*[cos(theta); sin(theta); 0];  %Position
V_pf = sqrt(mu/p)*[-sin(theta); (e+cos(theta)); 0];       %Velocity in 
                                                          %perifocal frame
R1p = [ cos(Omega)  -sin(Omega)  0;                %Calculate change of frame
         sin(Omega)  cos(Omega)  0;                %matrices with the 3
             0        0     1];                    %subsequent rotations
  
R2p = [1       0          0;
       0   cos(i)  -sin(i);
        0  sin(i)  cos(i)];


R3p = [ cos(w)  -sin(w)    0 ;
         sin(w)  cos(w)  0;
           0       0     1];
       
% Rotation matrix
T_pf2ge = R1p*R2p*R3p;                             %Obtain pf to fe matrix

r = T_pf2ge*R_pf;                        %Obtain the cartesian state vector
v = T_pf2ge*V_pf;                        %in geocentric frame
end

