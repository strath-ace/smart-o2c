function [FunVal] = SSFunctionFMINCON(y,thetaVec,trueThetaStep,TOF)

%SSFunctionFMINCON (Spherical Shaping Function)
% Function that calculates the DeltaTOF, the function to reduce to zero
% finding the root x (the vector of parameters of the spherical shaping)
% It is given as input to fmincon, with x as parameter to solve 
%
% INPUT
% y: parameters to find in the order (a0,a1,a2,a3,..,a6,b0,b3)
% thetaVec: vector with angles theta
% trueThetaStep: theta step for the numerical quadrature
% TOF: desired time of flight
%
%OUTPUT
% FunVal: value of the function DeltaTOF = abs(TOFcomp-TOF) to drive to zero
%         The absolute value is necessary due to the use in a minimization
%         algorithm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Niccolo' Gastaldello, October 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate the function of which to find the root

% Write down the state vector (with a2)
x(1:2,1) = y(1:2);
x(3:10,1) = y(4:11);
a2 = y(3);

% Add check of isDpos to avoid generation of imaginary numbers and
% subsequent generation of negative DV (do this taking only the real part
% of TOFcomp)

if (isDpos == 1)
    % Calculate T' and T''
    [TPrime1,TPrime2] = timeFunc(R,D,RPrime1,DPrime);
    
    % Calculate the time of flight with a Simpson quadrature
    TOFcomp = (trueThetaStep/6)*(TPrime1(1)+2*sum(TPrime1(3:2:end-1))+4*sum(TPrime1(2:2:end-1))+TPrime1(end));
    
    % Write the function to drive to zero and find the root (a2)
    FunVal = TOFcomp - TOF;
    
else
    FunVal = 0;
end

end

