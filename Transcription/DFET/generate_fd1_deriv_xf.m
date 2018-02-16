function M = generate_fd1_deriv_xf(x0,u0,t0,xf,uf,tf,x,u,t,static,scales,constants,fun,h)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%
% Auto generate forward finite difference approximation of the jacobian of 
% a given function wrt to final states xf. First order accuracy.

dx0 = fun(x0,u0,t0,xf,uf,tf,x,u,t,static,scales,constants);
dx0 = dx0(:,1);

M = zeros(length(dx0),length(x));

for jj = 1:length(x)
   
    x_temp = xf;
    x_temp(jj) = x_temp(jj)+h;
  
    dxtemp = fun(x0,u0,t0,x_temp,uf,tf,x,u,t,static,scales,constants);
    dxtemp = dxtemp(:,1);
  
    M(:,jj) = (dxtemp-dx0)/h;
    
end

end