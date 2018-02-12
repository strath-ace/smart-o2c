function [A, P, E, ERROR, VI, VF, TPAR, THETA] = lambertMR(RI, RF, TOF, MU, orbitType, Nrev, Ncase, optionsLMR)
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%-----------Copyright (C) 2016 University of Strathclyde-------------
%
%
%
%% lambertMR: Lambert's problem solver for all possible transfers:
%  
% 1- zero-revolution (for all possible types of orbits: circles, ellipses, 
%                  parabolas and hyperbolas)
% 2- multirevolution case
% 3- inversion of the motion
% 
% [A, P, E, ERROR, VI, VF, TPAR, THETA] = lambertMR(RI, RF, TOF, MU, orbitType, Nrev, Ncase)
% 
% 1- ZERO-REVOLUTION LAMBERT'S PROBLEM
% 
% For the solution of Lambert's problem with number of revolution = 0 the
% subroutine by Chris D'Souza is included here.
% This subroutine is a Lambert Algorithm which given two radius vectors and
% the time to get from one to the other, it finds the orbit connecting the
% two. It solves the problem using a new algorithm developed by R. Battin.
% It solves the Lambert problem for all possible types of orbits (circles, 
% ellipses, parabolas and hyperbolas).
% The only singularity is for the case of a transfer angle of 360 degrees, 
% which is a rather obscure case.
% It computes the velocity vectors corresponding to the given radius
% vectors except for the case when the transfer angle is 180 degrees
% in which case the orbit plane is ambiguous (an infinite number of
% transfer orbits exist).
% 
% 2- MULTIREVOLUTION LAMBERT'S PROBLEM
% 
% For the solution of Lambert's problem with Nrev>0 number of revolution, 
% Battin's formulation has been extended to accomodate N-revolution
% transfer orbits, by following the paper: "USING BATTIN METHOD TO OBTAIN 
% MULTIPLE-REVOLUTION LAMBERT'S SOLUTIONS" by Shen and Tsiotras.
% 
% When Nrev>0 the possible orbits are just ellipses.
% If 0<=Nrev<=Nmax, there are two Nrev-revolution transfer orbits.
% These two transfer orbits have different semi-major axis and they may be 
% all combinations of large-e and small-e transfer orbits.
% The Original Successive Substitution Method by Battin converges to one of
% the two possible solution with a viable initial guest, however it
% diverges from the other one. Then a Reversed Successive Substitution is
% used to converge to the second solution.
% A procedure is implemented in order to guarantee to provide initial
% guesses in the convergence region. If Nrew exceeds the maximum number of
% revolution an ERROR is given:
% warning('off', 'lambertMR:SuccessiveSubstitutionDiverged') to take out the
% warnings or use optionsLMR(1) = 0.
% 
% 
% 3- INVERSION OF THE MOTION
% 
% Direct or retrograde option can be selected for the transfer
% 
% The algorithm computes the semi-major axis, the parameter (semi-latus 
% rectum), the eccentricity and the velocity vectors.
% 
% NOTE: If ERROR occurs or the 360 or 180 degree transfer case is 
% encountered. 
% 
%% Inputs:
%      RI(3) = a three element array containing the initial position
%              vector [L]
%      RF(3) = a three element array containing the final position vector
%              [L]
%      TOF = transfer time, time of flight [T]
%      MU = planetary constant of the planet (mu = mass * G) [L^3/T^2]
%      orbitType = logical variable defining whether transfer is
%                  0: direct transfer from R1 to R2 (counterclockwise)
%                  1: retrograde transfer from R1 to R2 (clockwise)
%      Nrev = number of revolutions.
%             if Nrev = 0 ZERO-REVOLUTION transfer is calculated
%             if Nrev > 0 two transfers are possible. Ncase should be
%                         defined to select one of the two.
%      Ncase = logical variable defining the small-a or large-a option in 
%              case of Nrev>0:
%                  0: small-a option
%                  1: large-a option
%      optionsLMR = lambertMR options:
%                   optionsLMR(1) = display options:
%                                   0: no display
%                                   1: warnings are displayed only when
%                                      the algorithm does not converge
%                                   2: full warnings displayed
% 
%% Outputs:
%      A = Semi-major axis of the transfer orbit [L]
%      P = Semi-latus rectum of the transfer orbit [L]
%      E = Eccentricity of the transfer orbit
%      ERROR = Error flag
%              0:	No error
%              1:	Error, routine failed to converge
%              -1:	180 degrees transfer
%              2:  360 degrees transfer
%              3:  the algorithm doesn't converge because the number of
%                  revolutions is bigger than Nrevmax for that TOF
%              4:  Routine failed to converge, maximum number of
%                  iterations exceeded.
%      VI(3) = a three element array containing the initial velocity 
%              vector [L/T]
%      VT(3) = a three element array containing the final velocity vector
%              [L/T]
%      TPAR = parabolic flight time between RI and RF [T]
%      THETA = transfer angle [radians]
% 
% NOTE: The semi-major axis, positions, times, & gravitational parameter
%      must be in compatible units.
% 
% functions called: qck, cubicN (added at the bottom of this file)
% 
% issues to be solved: - 180 degrees transfer indetermination
%                     - 360 degrees transfer singularity
%                     - Nmax number of max revolution for a given TOF:
%                       work in progress - Camilla Colombo
% f
% evolution of: lambert.m (direct transfer, without multi-revolution
%                         option)
%                          PROGRAMMER:    Chris D'Souza
%                          DATE:          January 20, 1989
%                          VERIFIED BY:   Darrel Monroe, 10/25/90
%% Author(s):
% - Camilla Colombo - 10/11/2006
%                  - 13/11/2006 - added ERROR = 3 if Nrev > NrevMAX
%                  - 21/11/2006 - added another case of ERROR = 3 (index
%                                 N3) corresponding to the limit case when
%                                 small-a solution = large-a solution. No
%                                 solution is given in this case.
%                  - 06/08/2007 - optionsLMR added as an input
%                  - 28/11/2007 - minor changes
% 
% - Matteo Ceriotti - 29/01/2009
%  - Introduced variable for maximum number of iterations nitermax.
%  - Corrected final check on maximum number of iterations exceeded, from
%      "==" to ">=" (if N1 >= nitermax || N >= nitermax).
%  - Increased maxumum number of iterations to 2000, not to lose some
%      solutions.
%  - In OSS loop, added check for maximum number of iterations exceeded, 
%      which then sets checkNconvOSS = 0.
%  - Changed the way of coumputing X given Y1 in RSS. Now the
%      Newton-Raphson method with initial guess suggested by Shen, 
%      Tsiotras is used. This should guarantee convergence without the
%      need of an external zero finder (fsolve).
%  - Changed absolute tolerance into relative tolerance in all loops X0-X.
%      Now the condition is: while "abs(X0-X) >= abs(X)*TOL+TOL".
%  - Added return immediately when any error is detected.
%  - Moved check on 4*TOF*LAMBDA==0 after computing LAMBDA.
%  - Moved check on THETA==0 || THETA==2*PI after computing THETA.
%  - Added error code 4 (number of iterations exceeded).
%  - Removed variable Nwhile, as not strictly needed.
%  - Removed variable PIE=pi.
% 
% - Revised by Matteo Ceriotti - 29/01/2009.
% 
% Note: Please if you have got any changes that you would like to be done, 
%      do not change the function and let me know, because I am working on
%      it in order to improve it and solve the remaining issues. Camilla
% 
% ------------------------- - SpaceART Toolbox - --------------------------

% Check inputs
if nargin < 8
    optionsLMR = 0;
    if nargin < 6
        Nrev = 0;
        if nargin < 5
            orbitType = 0;
            if nargin < 4
                error('Not enough input arguments. See lambertMR.');
            end
        end
    end
end

nitermax = 2000; % Maximum number of iterations for loops
TOL = 1e-14;

TWOPI=2*pi;

% RESET
A=0;P=0;E=0;VI=[0, 0, 0];VF=[0, 0, 0];

% C
% C     ***  Compute the vector magnitudes and various
% C     ***  cross and dot products
% C

RIM2   = dot(RI, RI);
RIM    = sqrt(RIM2);
RFM2   = dot(RF, RF);
RFM    = sqrt(RFM2);
CTH    = dot(RI, RF)/(RIM*RFM);
CR     = cross(RI, RF);
STH    = norm(CR)/(RIM*RFM);

% *** choose angle for up angular momentum ***
switch orbitType
    case 0 % direct transfer
        if CR(3) < 0 
            STH = -STH;
        end
    case 1 % retrograde transfer
        if CR(3) > 0 
            STH = -STH;
        end
    otherwise
		error('% d is not an allowed orbitType', orbitType);
end
        
THETA  = qck(atan2(STH, CTH));

% if abs(THETA - pi) >= 0.01
if THETA == TWOPI || THETA==0
    ERROR = 2;
    A=0; P=0; E=0; VI=[0, 0, 0]; VF=[0, 0, 0]; TPAR=0; THETA=0;
    return
end

B1     = sign(STH); if STH == 0; B1 = 1; end;

% C
% C     ***  Compute the chord and the semi-perimeter
% C
C= sqrt(RIM2 + RFM2 - 2*RIM*RFM*CTH);
S= (RIM + RFM + C)/2;
BETA   = 2*asin(sqrt((S-C)/S));
PMIN   = TWOPI*sqrt(S^3/(8*MU));
TMIN   = PMIN*(pi-BETA+sin(BETA))/(TWOPI);
LAMBDA = B1*sqrt((S-C)/S);

if 4*TOF*LAMBDA == 0
    ERROR = -1;
    A=0; P=0; E=0; VI=[0, 0, 0]; VF=[0, 0, 0]; TPAR=0; THETA=0;
    return
end

% C
% C     ***  Compute L carefully for transfer angles less than 5 degrees
% C

if THETA*180/pi <= 5
   W   = atan((RFM/RIM)^.25) - pi/4;
   R1  = (sin(THETA/4))^2;
   S1  = (tan(2*W))^2;
   L   = (R1+S1)/(R1+S1+cos(THETA/2));
else
   L   = ((1-LAMBDA)/(1+LAMBDA))^2;
end

M= 8*MU*TOF^2/(S^3*(1+LAMBDA)^6);
TPAR   = (sqrt(2/MU)/3)*(S^1.5-B1*(S-C)^1.5);
L1     = (1 - L)/2;

CHECKFEAS = 0;
N1 = 0;
N = 0;

if Nrev == 0
    % C
    % C     ***  Initialise values of y, n, and x
    % C

    Y= 1;
    N= 0;
    N1=0;
    ERROR  = 0;
    % CHECKFEAS=0;

    if (TOF-TPAR) <= 1e-3
        X0  = 0;
    else
        X0  = L;
    end

    X= -1.e8;

    % C
    % C     ***  Begin iteration
    % C
    
    % Matteo Ceriotti, 26/01/2009
    % Changed absolute tolerance into relative tolerance here below.
    while (abs(X0-X) >= abs(X)*TOL+TOL) && (N <= nitermax)
        N   = N+1;
        X   = X0;
        ETA = X/(sqrt(1+X) + 1)^2;
        CHECKFEAS=1;

        % C
        % C  ***  Compute x by means of an algorithm devised by
        % C  ***  Gauticci for evaluating continued fractions by the
        % C  ***  'Top Down' method
        % C

        DELTA = 1;
        U     = 1;
        SIGMA = 1;
        M1    = 0;

        while abs(U) > TOL && M1 <= nitermax
            M1    = M1+1;
            GAMMA = (M1 + 3)^2/(4*(M1+3)^2 - 1);
            DELTA = 1/(1 + GAMMA*ETA*DELTA);
            U     = U*(DELTA - 1);
            SIGMA = SIGMA + U;
        end

        C1 = 8*(sqrt(1+X)+1)/(3+1/(5 + ETA + (9*ETA/7)*SIGMA));

        % C
        % C  ***  Compute H1 and H2
        % C

        if N == 1
            DENOM = (1 + 2*X + L)*(3*C1 + X*C1 +4*X);
            H1 = (L+X)^2*(C1 + 1 + 3*X)/DENOM;
            H2 = M*(C1+X-L)/DENOM;
        else
            QR = sqrt(L1^2 + M/Y^2);
            XPLL = QR - L1;
            LP2XP1 = 2*QR;
            DENOM = LP2XP1*(3*C1 + X*C1+4*X);
            H1 = ((XPLL^2)*(C1 + 1 + 3*X))/DENOM;
            H2 = M*(C1+X-L)/DENOM;
        end
        
        B = 27*H2/(4*(1+H1)^3);
        U = -B/(2*(sqrt(B+1)+1));

        % C
        % C  ***  Compute the continued fraction expansion K(u)
        % C  ***  by means of the 'Top Down' method
        % C

        % Y can be computed finding the roots of the formula and selecting
        % the real one:
        % y^3 - (1+h1)*y^2 - h2 = 0     (7.113) Battin
        % 
        % Ycami_ = roots([1 -1-H1 0 -H2])
        % kcami = find( abs(imag(Ycami_)) < eps );
        % Ycami = Ycami_(kcami)

        DELTA = 1;
        U0 = 1;
        SIGMA = 1;
        N1 = 0;

        while N1 < nitermax && abs(U0) >= TOL
            if N1 == 0
                GAMMA = 4/27;
                DELTA = 1/(1-GAMMA*U*DELTA);
                U0 = U0*(DELTA - 1);
                SIGMA = SIGMA + U0;
            else
                for I8 = 1:2
                    if I8 == 1
                        GAMMA = 2*(3*N1+1)*(6*N1-1)/(9*(4*N1 - 1)*(4*N1+1));
                    else
                        GAMMA = 2*(3*N1+2)*(6*N1+1)/(9*(4*N1 + 1)*(4*N1+3));
                    end
                    DELTA = 1/(1-GAMMA*U*DELTA);
                    U0 = U0*(DELTA-1);
                    SIGMA = SIGMA + U0;
                end
            end

            N1 = N1 + 1;
        end

        KU = (SIGMA/3)^2;
        Y = ((1+H1)/3)*(2+sqrt(B+1)/(1-2*U*KU));    % Y = Ycami
        
        X0 = sqrt(((1-L)/2)^2+M/Y^2)-(1+L)/2;
        % fprintf('n= % d, x0=% .14f\n', N, X0);
    end
    
% MULTIREVOLUTION
elseif (Nrev > 0) && (4*TOF*LAMBDA~=0) % (abs(THETA)-pi > 0.5*pi/180)

    checkNconvRSS = 1;
    checkNconvOSS = 1;
    N3 = 1;
    
    while N3 < 3
        
        if Ncase == 0 || checkNconvRSS == 0

            % - Original Successive Substitution -
            % always converges to xL - small a

            % C
            % C     ***  Initialise values of y, n, and x
            % C

            Y= 1;
            N= 0;
            N1=0;
            ERROR = 0;
            % CHECKFEAS = 0;
%            if (TOF-TPAR) <= 1e-3
%                X0 = 0;
%            else
            if checkNconvOSS == 0
                X0 = 2*X0;
                checkNconvOSS = 1;
                % see p. 11 USING BATTIN METHOD TO OBTAIN 
                % MULTIPLE-REVOLUTION LAMBERT'S SOLUTIONS - Shen, Tsiotras
            elseif checkNconvRSS == 0;
                % X0 is taken from the RSS
            else
                X0 = L;
            end

            X = -1.e8;

            % C
            % C     ***  Begin iteration
            % C
			
            % Matteo Ceriotti, 26/01/2009
            % Changed absolute tolerance into relative tolerance here
            % below.
            while (abs(X0-X) >= abs(X)*TOL+TOL) && (N <= nitermax)
                N   = N+1;
                X   = X0;
                ETA = X/(sqrt(1+X) + 1)^2;
                CHECKFEAS = 1;

                % C
                % C  ***  Compute x by means of an algorithm devised by
                % C  ***  Gauticci for evaluating continued fractions by the
                % C  ***  'Top Down' method
                % C

                DELTA = 1;
                U     = 1;
                SIGMA = 1;
                M1    = 0;

                while abs(U) > TOL && M1 <= nitermax
                    M1    = M1+1;
                    GAMMA = (M1 + 3)^2/(4*(M1+3)^2 - 1);
                    DELTA = 1/(1 + GAMMA*ETA*DELTA);
                    U     = U*(DELTA - 1);
                    SIGMA = SIGMA + U;
                end

                C1 = 8*(sqrt(1+X)+1)/(3+1/(5 + ETA + (9*ETA/7)*SIGMA));

                % C
                % C  ***  Compute H1 and H2
                % C

                if N == 1
                    DENOM = (1 + 2*X + L)*(3*C1 + X*C1 +4*X);
                    H1 = (L+X)^2*(C1 + 1 + 3*X)/DENOM;
                    H2 = M*(C1+X-L)/DENOM;
                else
                    QR = sqrt(L1^2 + M/Y^2);
                    XPLL = QR - L1;
                    LP2XP1 = 2*QR;
                    DENOM = LP2XP1*(3*C1 + X*C1+4*X);
                    H1 = ((XPLL^2)*(C1 + 1 + 3*X))/DENOM;
                    H2 = M*(C1+X-L)/DENOM;
                end

                H3 = M*Nrev*pi/(4*X*sqrt(X));
                H2 = H3+H2;

                B = 27*H2/(4*(1+H1)^3);
                U = -B/(2*(sqrt(B+1)+1));

                % C
                % C  ***  Compute the continued fraction expansion K(u)
                % C  ***  by means of the 'Top Down' method
                % C

                % Y can be computed finding the roots of the formula and selecting
                % the real one:
                % y^3 - (1+h1)*y^2 - h2 = 0     (7.113) Battin
                % 
                % Ycami_ = roots([1 -1-H1 0 -H2])
                % kcami = find( abs(imag(Ycami_)) < eps );
                % Ycami = Ycami_(kcami)

                DELTA = 1;
                U0 = 1;
                SIGMA = 1;
                N1 = 0;

                while N1 < nitermax && abs(U0) >= TOL
                    if N1 == 0
                        GAMMA = 4/27;
                        DELTA = 1/(1-GAMMA*U*DELTA);
                        U0 = U0*(DELTA - 1);
                        SIGMA = SIGMA + U0;
                    else
                        for I8 = 1:2
                            if I8 == 1
                                GAMMA = 2*(3*N1+1)*(6*N1-1)/(9*(4*N1 - 1)*(4*N1+1));
                            else
                                GAMMA = 2*(3*N1+2)*(6*N1+1)/(9*(4*N1 + 1)*(4*N1+3));
                            end
                            DELTA = 1/(1-GAMMA*U*DELTA);
                            U0 = U0*(DELTA-1);
                            SIGMA = SIGMA + U0;
                        end
                    end

                    N1 = N1 + 1;
                end

                KU = (SIGMA/3)^2;
                Y = ((1+H1)/3)*(2+sqrt(B+1)/(1-2*U*KU));	% Y = Ycami
                if Y > sqrt(M/L)
                    if optionsLMR(1) == 2
                        warning('lambertMR:SuccessiveSubstitutionDiverged', ...
                                ['Original Successive Substitution is diverging\n'...
                                '-> Reverse Successive Substitution used to find the proper XO.\n']);
                    end
                    checkNconvOSS = 0;
                    break
                end
                
                X0 = sqrt(((1-L)/2)^2+M/Y^2)-(1+L)/2;
                % fprintf('N: % d X0: % .14f\n', N, X0);
            end
            
            % When 2 solutions exist (small and big a), the previous loop
            % must either converge or diverge because Y > sqrt(M/L) at some
            % point. Thus, the upper bound on the number of iterations
            % should not be necessary. Though, nothing can be said in the
            % case tof<tofmin and so no solution exist. In this case, an
            % upper bound on number of iterations could be needed.
            
            if N >= nitermax % Checks if previous loop ended due to maximum number of iterations
                if optionsLMR(1) == 2
                    warning('lambertMR:SuccessiveSubstitutionExceedMaxIter', ...
                            ['Original Successive Substitution exceeded max number of iteration\n'...
                            '-> Reverse Successive Substitution used to find the proper XO.\n']);
                end
                checkNconvOSS = 0;
            end
        end
        if (Ncase == 1 || checkNconvOSS == 0) && ~(checkNconvRSS == 0 && checkNconvOSS == 0)

            % - Reverse Successive Substitution -
            % always converges to xR - large a

            % C
            % C     ***  Initialise values of y, n, and x
            % C

            N = 0;
            N1 = 0;
            ERROR  = 0;
            % CHECKFEAS=0;
            if checkNconvRSS == 0;
                X0 = X0/2; % XL/2
                checkNconvRSS = 1;
                % see p. 11 USING BATTIN METHOD TO OBTAIN 
                % MULTIPLE-REVOLUTION LAMBERT'S SOLUTIONS - Shen, Tsiotras
            elseif checkNconvOSS == 0
                % X0 is taken from the OSS
            else
                X0 = L;
            end

            X = -1.e8;

            % C
            % C     ***  Begin iteration
            % C
            
            % Modified by: Matteo Ceriotti, 26/01/2009
            % Changed absolute tolerance into relative tolerance here
            % below.
            while (abs(X0-X) >= abs(X)*TOL+TOL) && (N <= nitermax)
                N = N+1;
                X = X0;
                CHECKFEAS=1;

                Y = sqrt(M/((L+X)*(1+X))); % y1 in eq. (8a) in Shen, Tsiotras

                if Y < 1
                    if optionsLMR(1) == 2
                        warning('lambertMR:SuccessiveSubstitutionDiverged', ...
                                ['Reverse Successive Substitution is diverging\n' ...
                                '-> Original Successive Substitution used to find the proper XO.\n']);
                    end
                    checkNconvRSS = 0;
                    break
                end
                
                % Modified by: Matteo Ceriotti, 27/01/2009
                % This is the Newton-Raphson method suggested by USING
                % BATTIN METHOD TO OBTAIN MULTIPLE-REVOLUTION LAMBERT'S
                % SOLUTIONS - Shen, Tsiotras
                
                % To assure the Newton-Raphson method to be convergent
                Erss = 2*atan(sqrt(X));
                while h_E(Erss, Y, M, Nrev) < 0
                    Erss = Erss/2;
                end
                
                Nnew = 1;
                Erss_old = -1.e8;
                
                % The following Newton-Raphson method should always
                % converge, given the previous first guess choice, 
                % according to the paper. Therefore, the condition on
                % number of iterations should not be neccesary. It could be
                % necessary for the case tof < tofmin.
                while (abs(Erss-Erss_old) >= abs(Erss)*TOL+TOL) && Nnew < nitermax
                    Nnew = Nnew+1;
                    [h, dh] = h_E(Erss, Y, M, Nrev);
                    Erss_old = Erss;
                    Erss = Erss - h/dh;
                    % fprintf('Nnew: % d Erss: % .16f h_E: % .16f\n', Nnew, Erss, h);
                end
                if Nnew >= nitermax
                    if optionsLMR(1) ~= 0
                        warning('lambertMR:NewtonRaphsonIterExceeded', 'Newton-Raphson exceeded max iterations.\n');
                    end
                end
                X0 = tan(Erss/2)^2;
            end
        end
        if checkNconvOSS == 1 && checkNconvRSS == 1
            break
        end
        
        if checkNconvRSS == 0 && checkNconvOSS == 0
            if optionsLMR ~=0
                warning('lambertMR:SuccessiveSubstitutionDiverged', ...
                        ['Both Original Successive Substitution and Reverse ' ...
                        'Successive Substitution diverge because Nrev > NrevMAX.\n' ...
                        'Work in progress to calculate NrevMAX.\n']);
            end
            ERROR = 3;
            A=0; P=0; E=0; VI=[0, 0, 0]; VF=[0, 0, 0]; TPAR=0; THETA=0;
            return
        end
        
        N3 = N3+1;
    end
    
    if N3 == 3
        if optionsLMR ~=0
            warning('lambertMR:SuccessiveSubstitutionDiverged', ...
                    ['Either Original Successive Substitution or Reverse ' ...
                    'Successive Substitution is always diverging\n' ...
                    'because Nrev > NrevMAX or because large-a solution = small-a solution (limit case).\n' ...
                    'Work in progress to calculate NrevMAX.\n']);
        end
        ERROR = 3;
        A=0; P=0; E=0; VI=[0, 0, 0]; VF=[0, 0, 0]; TPAR=0; THETA=0;
        return
    end
end

% C
% C     ***  Compute the velocity vectors
% C

if CHECKFEAS == 0
    ERROR = 1;
    A=0; P=0; E=0; VI=[0, 0, 0]; VF=[0, 0, 0]; TPAR=0; THETA=0;
    return
end

if N1 >= nitermax || N >= nitermax
    ERROR = 4;
    if optionsLMR ~=0
        disp('Lambert algorithm has not converged, maximum number of iterations exceeded.');
    end
    A=0; P=0; E=0; VI=[0, 0, 0]; VF=[0, 0, 0]; TPAR=0; THETA=0;
    return
end

CONST = M*S*(1+LAMBDA)^2;
A = CONST/(8*X0*Y^2);

R11 = (1 + LAMBDA)^2/(4*TOF*LAMBDA);
S11 = Y*(1 + X0);
T11 = (M*S*(1+LAMBDA)^2)/S11;

VI(1:3) = -R11*(S11*(RI(1:3)-RF(1:3))-T11*RI(1:3)/RIM);
VF(1:3) = -R11*(S11*(RI(1:3)-RF(1:3))+T11*RF(1:3)/RFM);

P = (2*RIM*RFM*Y^2*(1+X0)^2*sin(THETA/2)^2)/CONST;
E = sqrt(1 - P/A);

return

% -------------------------------------------------------------------------



function [angle]=qck(angle)

%% ***********************************************************************
% *     This function takes any angle and reduces it, if necessary, 
% *     so that it lies in the range from 0 to 2 PI radians.
% *
% *     INPUTS TO THE FUNCTION
% *          ANGLE   = The ange to be reduced (in radians)
% *
% *     OUTPUTS FROM THE FUNCTION
% *          QCK     = The angle reduced, if necessary, to the range
% *                    from 0 to 2 PI radians (in radians)
% *
% *     MISSION PLANNING SUBROUTINES AND FUNCTIONS CALLED
% *          PI
% *
% *     PROGRAMMER:    W.T. Fowler
% *
% *     DATE:          July, 1978
% *
% *     VERIFIED BY:   Darrel Monroe, 8/20/90
% *
% ************************************************************************

twopi = 2*pi;
 
diff = twopi * (fix(angle/twopi) + min([0, sign(angle)]));

angle = angle -diff;

return

function [h, dh] = h_E(E, y, m, Nrev)
% Equation of multirevolution Lambert's problem h = h(E).
% 
%  h(E) = (Nrev*pi + E - sin(E)) / tan(E/2)^3 - 4/m * (y^3 - y^2)
% 
% See: "USING BATTIN METHOD TO OBTAIN MULTIPLE-REVOLUTION LAMBERT'S 
%     SOLUTIONS", Shen, Tsiotras, pag. 12
% 
% INPUT
%  E, y, m, Nrev   See paper for detailed description.
% 
% OUTPUT
%  h               Value of h(E).
%  dh              Value of dh(E)/dE.
% 
% Matteo Ceriotti, 27/01/2009

tanE2 = tan(E/2);
h = (Nrev*pi + E - sin(E)) / tanE2^3 - 4/m * (y^3 - y^2);

if nargout > 1  % two output arguments
    % h'(E)
    dh = (1-cos(E))/tanE2^3 - 3/2*(Nrev*pi+E-sin(E))*sec(E/2)^2 / tanE2^4;
end

return
