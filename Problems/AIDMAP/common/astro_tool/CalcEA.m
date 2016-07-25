% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%-----------Copyright (C) 2016 University of Strathclyde-------------
%


% Solves for eccentric anomaly, E from a given mean anomaly, M
% and eccentricty, e.  Performs a simple Newton-Raphson iteration
% 
% M must be in radians and E is returned in radians.
% 
% Default tolerances is 10^-8 rad.
% 
% E = CalcEA(M, e) uses default tolerances
% 
% E = CalcEA(M, e, tol) will use a user specified tolerance, tol
% 
% Richard Rieber
% 1/23/2005
% E-mail problems/suggestions rrieber@gmail.com

function E = CalcEA(M, e)

tol = 10^-8;

Etemp = M;
ratio = 1;

while abs(ratio) > tol
    
    f_E = Etemp - e*sin(Etemp) - M;
    
    f_Eprime = 1 - e*cos(Etemp);
    
    ratio = f_E/f_Eprime;
    
    if abs(ratio) > tol
        Etemp = Etemp - ratio;
    else
        E = Etemp;
    end
end