function dx = Two_burns_modified_state_equations(x,u,time,static,scales,constants)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%

% states:   x(1) = p
%           x(2) = f
%           x(3) = g
%           x(4) = h
%           x(5) = k
%           x(6) = L
%           x(7) = m

%
% controls: u(1) = theta
%           u(2) = psi

%% constants

xscale = scales.xscale;
uscale = scales.uscale;
tscale = scales.tscale;
staticscale = scales.other_vars_scale;

%% de normalising input

x = x.*xscale;
u = u.*uscale;                  
time = time.*tscale;
static = static.*staticscale;

%% the actual dynamics

p = x(1);
f = x(2);
g = x(3);
h = x(4);
k = x(5);
L = x(6)*pi/180;

q = 1+f*cos(L)+g*sin(L);

pm = (p/constants.mu)^0.5;
r = p/q;
a2 = h^2-k^2;
s2 = 1+(h^2+k^2);

A = [ 0 2*p/q*pm 0;
      pm*sin(L) pm/q*((q+1)*cos(L)+f) -pm*g/q*(h*sin(L)-k*cos(L));
     -pm*cos(L) pm/q*((q+1)*sin(L)+g)  pm*f/q*(h*sin(L)-k*cos(L));
     0 0 pm*s2*cos(L)/(2*q);
     0 0 pm*s2*sin(L)/(2*q);
     0 0 pm/q*(h*sin(L)-k*cos(L))];
 
b = [0;
     0;
     0;
     0;
     0;
     (constants.mu*p)^0.5*(q/p)^2];

rvett = [r/s2*(cos(L)+a2*cos(L)+2*h*k*sin(L));
         r/s2*(sin(L)-a2*sin(L)+2*h*k*cos(L));
         2*r/s2*(h*sin(L)-k*cos(L))];
 
vvett = [-1/(s2*pm) * (   sin(L) + a2*sin(L) - 2*h*k*cos(L) + g - 2*f*h*k + a2*g);
         -1/(s2*pm) * (  -cos(L) + a2*cos(L) + 2*h*k*sin(L) - f + 2*g*h*k + a2*f);
         2/(s2*pm) *  ( h*cos(L) + k*sin(L)  + f*h + g*k)];

rxv = cross(rvett,vvett);
vxr = cross(vvett,rvett);

Qr = [rvett/norm(rvett) cross(rxv,rvett)/norm(cross(rxv,rvett)) rxv/norm(rxv)];
Qv = [vvett/norm(vvett) vxr/norm(vxr) cross(vvett/norm(vvett),vxr/norm(vxr))];

T = constants.Tc*u(3)*constants.g0/x(7)*Qv*[cos(u(1)*pi/180)*cos(u(2)*pi/180); cos(u(1)*pi/180)*sin(u(2)*pi/180); sin(u(1)*pi/180)];

delta = Qr'*T;
     
dx = [A*delta+b; -constants.Tc*u(3)/constants.Isp];

%% renormalise output

dx = dx*tscale./(xscale);

end