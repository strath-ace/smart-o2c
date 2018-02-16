function [cart] = equinoctial_to_cartesian(mod,mu)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%

p = mod(1);
f = mod(2);
g = mod(3);
h = mod(4);
k = mod(5);
L = mod(6);

q = 1+f*cos(L)+g*sin(L);
r = p/q;
s2 = 1+(h^2+k^2);
a2 = h^2-k^2;
pm = (p/mu)^0.5;

rvett = [r/s2*(cos(L)+a2*cos(L)+2*h*k*sin(L));
         r/s2*(sin(L)-a2*sin(L)+2*h*k*cos(L));
         2*r/s2*(h*sin(L)-k*cos(L))];
 
vvett = [-1/(s2*pm) * (   sin(L) + a2*sin(L) - 2*h*k*cos(L) + g - 2*f*h*k + a2*g);
         -1/(s2*pm) * (  -cos(L) + a2*cos(L) + 2*h*k*sin(L) - f + 2*g*h*k + a2*f);
         2/(s2*pm) *  ( h*cos(L) + k*sin(L)  + f*h + g*k)];
     
cart = [rvett' vvett'];

end