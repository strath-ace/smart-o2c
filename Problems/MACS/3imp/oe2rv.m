% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function rv = oe2rv(orb_elem,planet);

a = orb_elem(1);                                                            % semimajor axis
e = orb_elem(2);                                                            % eccentricity
i = orb_elem(3);                                                            % inclination of orbital plan
th = orb_elem(4);                                                           % true anomaly
W = orb_elem(5);                                                            % right ascension of the ascending node
w = orb_elem(6);                                                            % perigee argument

mu = planet(1);

f = th-w;
p = a*(1-e^2);
h = sqrt(p*mu);

r = p/(1+e*cos(f));

R = r.*[cos(W)*cos(th)-sin(W)*sin(th)*cos(i);
        sin(W)*cos(th)+cos(W)*sin(th)*cos(i);
        sin(th)*sin(i)];
    
V = (mu/h).*[-(cos(W)*(sin(th)+e*sin(w))+sin(W)*(cos(th)+e*cos(w))*cos(i));
             -(sin(W)*(sin(th)+e*sin(w))-cos(W)*(cos(th)+e*cos(w))*cos(i));
             (cos(th)+e*cos(w))*sin(i)];
rv = [R;V];
