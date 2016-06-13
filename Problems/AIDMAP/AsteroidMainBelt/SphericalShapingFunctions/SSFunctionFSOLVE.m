function [FunVal] = SSFunctionFSOLVE(a2,inCond,A,B,thetaVec,trueThetaStep,TOF)

%SSFunctionFSOLVE (Spherical Shaping Function)
% Function that calculates the DeltaTOF, the function to reduce to zero
% finding the root a2 (the 3rd parameter of the spherical shaping)
% It is given as input to fsolve, with only a2 as parameter to solve (the
% other are found exactly from the resolution of the linear system)
%
% INPUT
% a2: parameter to find
% inCond: struct with initial conditions
% A,B: matrix and vector of the linear system
% thetaVec: vector with angles theta
% trueThetaStep: theta step for the numerical quadrature
% TOF: desired time of flight
%
%OUTPUT
% FunVal: value of the function DeltaTOF = TOFcomp-TOF to drive to zero

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Niccolo' Gastaldello, October 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate the function of which to find the root

% Write down the initial conditions
Theta_i = inCond.ti;
Theta_f = inCond.tf;
R_i = inCond.ri;
R_f = inCond.rf;

% Construct the vector Aa2
Aa2 = [a2*Theta_i^2; a2*Theta_f^2; 0; 0; a2*2*Theta_i; a2*2*Theta_f; 0; 0; -a2*2*R_i^2; -a2*2*R_f^2];

% Solve for the parameters x
x = A\(B-Aa2);

% Calculate the functions and their derivatives
[R,Phi,RPrime1,PhiPrime1,RPrime2,PhiPrime2,RPrime3,PhiPrime3,D,DPrime,isDpos] = baseFuncDeriv(thetaVec,[x;a2]);

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

