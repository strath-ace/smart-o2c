function y = CloneIC(x_0,structure,t)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%
% Clones initial condition over collocation nodes, with linear 
% extrapolation between the boundary conditions if both are given

num_times = length(t)-1*(structure.DFET==1)*(structure.state_order==0);
num_free_end_states = sum(structure.free_final_states);

if size(x_0,1)>1
   
    x_0 = x_0';
    
end

y = zeros(num_times+(num_free_end_states>0)*(structure.DFET==1),structure.num_eqs);
y(1,:) = x_0;

x_f = structure.x_f;

if size(x_f)~=size(x_0)
   
    x_f = x_f';
    
end

x_f(~structure.imposed_final_states) = x_0(~structure.imposed_final_states);

dxdt = (x_f-x_0)/(t(end)-t(1));

for i =2:num_times    

    %y(i,:) = x_0;
    y(i,:) = x_0+dxdt*(t(i)-t(1));
    
end

%this has to be improved/specialized
if structure.DFET==1 && (num_free_end_states>0) 

    y(end,:) = y(end-1,:);
    
end

end