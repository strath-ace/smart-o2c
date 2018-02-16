function [kep] = equinoctial_to_keplerian(mod)

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

e = (f^2+g^2)^0.5;
a = p/(1-e^2);
O = atan2(h,k);
i = 2*atan2(cos(O),h);
if e < eps

    w = 0;
    
else
    
    w = acos(f/e)-O;
    
end

ni = L-O-w;

kep = [a,e,i,O,w,ni];

end