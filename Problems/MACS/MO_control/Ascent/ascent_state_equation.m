% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%-----------Copyright (C) 2016 University of Strathclyde-------------
%


function M = ascent_state_equation(x,u,t)

g = 0.0016; 
T = 4e-3;

M = [x(2); 
    T*cos((u(1))); 
    x(4); 
    T*sin((u(1)))-g]; 

end