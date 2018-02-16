function [x0,u0,t0,xf,uf,tf,static] = get_phase_boundary_vals (x_in,structure)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%
% This function retrieves the NORMALISED boundary states x0 and xf, 
% controls u0 and uf, times t0 and tf, and the static variables.

[x0,xb,t0,tf,static] = get_boundary_and_static_vars(x_in,structure); % beware, xb here are the IMPOSED final states. ALL STATES ARE NORMALISED

% Removing static variables

x = x_in(~structure.static_vars);    
[~,~,xf] = extract_solution(x,structure,xb); % the actual NORMALISED values of xf

% mapping mixed states and controls of input x vector into more easily
% managed states and controls

loc_state = zeros((structure.state_order+1)*structure.num_eqs,structure.num_elems); % preallocation, will be needed to more easily map between the input vector (mixed states and controls) and the states
loc_u = zeros((structure.control_order+1)*structure.num_controls,structure.num_elems); % preallocation, will be needed to more easily map between the input vector (mixed states and controls) and the controls

startx = 1;

for i=1:structure.num_elems
    
    endx = startx+(structure.state_order+1)*structure.num_eqs-1;
    loc_state(:,i) = x(startx:endx);
    
    startx = endx+1*(structure.num_controls>0);
    endx = startx+(structure.control_order+1)*structure.num_controls-1*(structure.num_controls>0);
    
    loc_u(:,i) = x(startx:endx);
    
    startx = endx+1;
    
end

u0 = structure.control_valsl*loc_u(:,1);
uf = structure.control_valsr*loc_u(:,end);

end