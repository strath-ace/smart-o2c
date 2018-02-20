function g=g_fun(f,lambda,z,zstar)
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%-----------Copyright (C) 2016 University of Strathclyde-------------
%
%
%
%g=max(lambda.*abs(f-z));

% in case zstar==z we have a Nan... to avoid it, i don't normalise in that
% case (also because in the single objective case zstar always equals z)

f2 = (f-z)./(((zstar-z).*((zstar-z)~=0))+((zstar-z)==0));

[~,index] = max(lambda.*abs(f2));
g = lambda(index).*f2(index);

% assuming ||lambda|| = 1
%d1 = abs((f-z)*lambda');
%d2 = ((f-(z+d1*lambda))*(f-(z+d1*lambda))')^0.5;

%g = d1+5*d2;

return
