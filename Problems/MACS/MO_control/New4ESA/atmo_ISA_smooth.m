function [p,rho,c,T] = atmo_ISA_smooth(hm)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------

%function [p, atmo.dens, atmo.speed, atmo.temp] = atmo_ISA_smooth(hm)
%{
Smoothed semi-continuous function of the ISA atmospheric model
 atmo = atmo_ISA_smooth(h)

INPUT 
h = altitude above mean sea level (m)

OUTPUT
atmo = structure containing [pres, temp, density, sound] in SI units

Note: based on curve fitting, pressure is continuous, temperature continuous up to 86 km
    after which a static minimum value is imposed. 

(c) 2016, cFASTT
%}

h = hm.*1e-3; %fitting based on km

% pressure based on Gaussian
b = [2618.07756933778,2.71920591908197,6.84221679629987,485470.873628963,-32.1750294095617,22.4313026971687,258.132670001550,8.57839616821948,3.14058778995278,4.55750594929836e+15,-432.227390750371,85.5349847929330];
pres = b(1)*exp(-((h-b(2))/b(3))^2) + b(4)*exp(-((h-b(5))/b(6))^2) + ...
              b(7)*exp(-((h-b(8))/b(9))^2) + b(10)*exp(-((h-b(11))/b(12))^2);
p = max(0, pres);

% temperature based on fourier series
a = [5607.98771649688,11554.9833043799,-10053.3072049395,-8203.10684485628,-18188.0803970065,-18149.7018348519,5724.18950065030,3377.72970310996,13421.9260350188,7525.62142206319,-1599.91159138216,-572.123833720502,-3116.96671034414,-869.211521590768,135.248232816609,17.3109846416532,127.432888173767,0.0349697229334475];
w = a(18); 
T =   a(1) + a(2)*cos(h*w) + a(3)*sin(h*w) + ...
               a(4)*cos(2*h*w) + a(5)*sin(2*h*w) + a(6)*cos(3*h*w) + a(7)*sin(3*h*w) + ...
               a(8)*cos(4*h*w) + a(9)*sin(4*h*w) + a(10)*cos(5*h*w) + a(11)*sin(5*h*w) + ...
               a(12)*cos(6*h*w) + a(13)*sin(6*h*w) + a(14)*cos(7*h*w) + a(15)*sin(7*h*w) + ...
               a(16)*cos(8*h*w) + a(17)*sin(8*h*w);
T = max(186.87, T);
           
rho=p/(287.058*T);
c=sqrt(1.4*287.058*T);
           
return
           
%{
General model Fourier8: TEMPERATURE
     f(x) = 
               a0 + a1*cos(x*w) + b1*sin(x*w) + 
               a2*cos(2*x*w) + b2*sin(2*x*w) + a3*cos(3*x*w) + b3*sin(3*x*w) + 
               a4*cos(4*x*w) + b4*sin(4*x*w) + a5*cos(5*x*w) + b5*sin(5*x*w) + 
               a6*cos(6*x*w) + b6*sin(6*x*w) + a7*cos(7*x*w) + b7*sin(7*x*w) + 
               a8*cos(8*x*w) + b8*sin(8*x*w)
Coefficients (with 95% confidence bounds):
       a0 =        5608  (-1.687e+04, 2.808e+04)
       a1 =   1.155e+04  (-1.113e+04, 3.424e+04)
       b1 =  -1.005e+04  (-5.005e+04, 2.994e+04)
       a2 =       -8203  (-3.642e+04, 2.001e+04)
       b2 =  -1.819e+04  (-5.335e+04, 1.698e+04)
       a3 =  -1.815e+04  (-5.232e+04, 1.602e+04)
       b3 =        5724  (-1.053e+04, 2.198e+04)
       a4 =        3378  (-5484, 1.224e+04)
       b4 =   1.342e+04  (-1.078e+04, 3.762e+04)
       a5 =        7526  (-5268, 2.032e+04)
       b5 =       -1600  (-6769, 3569)
       a6 =      -572.1  (-3167, 2022)
       b6 =       -3117  (-8009, 1775)
       a7 =      -869.2  (-2102, 363.3)
       b7 =       135.2  (-743, 1014)
       a8 =       17.31  (-135.7, 170.3)
       b8 =       127.4  (-27.55, 282.4)
       w =     0.03497  (0.0314, 0.03854)

Goodness of fit:
  SSE: 478.7
  R-square: 0.9992
  Adjusted R-square: 0.9992
  RMSE: 0.6982



-------------
General model Gauss4: PRESSURE
     f(x) = 
              a1*exp(-((x-b1)/c1)^2) + a2*exp(-((x-b2)/c2)^2) + 
              a3*exp(-((x-b3)/c3)^2) + a4*exp(-((x-b4)/c4)^2)
Coefficients (with 95% confidence bounds):
       a1 =        2618  (2403, 2834)
       b1 =       2.719  (2.412, 3.026)
       c1 =       6.842  (6.692, 6.993)
       a2 =   4.855e+05  (1.291e+05, 8.419e+05)
       b2 =      -32.18  (-38.45, -25.9)
       c2 =       22.43  (20.77, 24.09)
       a3 =       258.1  (234.8, 281.5)
       b3 =       8.578  (8.525, 8.631)
       c3 =       3.141  (3.04, 3.241)
       a4 =   4.558e+15  (-1.457e+17, 1.548e+17)
       b4 =      -432.2  (-930.2, 65.75)
       c4 =       85.53  (41.56, 129.5)

Goodness of fit:
  SSE: 9740
  R-square: 1
  Adjusted R-square: 1
  RMSE: 3.14

%}


