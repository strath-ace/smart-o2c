function [x_0,x_f,t_0,t_f,static,ids] = get_boundary_and_static_vars(x_in,structure)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%
% This function retrieves the NORMALISED boundary states x0 and xf,
% times t0 and tf, the static variables, and the ids of free initial states

static = x_in(structure.other_vars);

if structure.imposed_t0
    
    t_0 = structure.t_0_norm;
    
else
    
    t_0 = x_in(structure.t0_vars);
    
end

if structure.imposed_tf
    
    t_f = structure.t_f_norm;
    
else
    
    t_f = x_in(structure.tf_vars);
    
end

ids = 1:length(structure.x_0);
ids = ids(~structure.imposed_initial_states);   % get which x_0 are coded in solution vector
x_0 = structure.x_0_norm;
x_0(ids) = x_in(structure.x0_vars);
x_f = structure.x_f_norm;

end