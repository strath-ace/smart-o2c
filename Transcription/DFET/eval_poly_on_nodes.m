function [state_eval_nodes,control_eval_nodes,t_nodes] = eval_poly_on_nodes(num_elems,state_basis,control_basis,t_nodes)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%
% PROBABLY OBSOLETE
% Evaluate polynomials of states and controls on the same given nodes

%enodes = sort(unique([[-1 state_nodes 1],control_nodes])); % including -1 and 1 artificially can generate problems afterwards

state_eval_nodes = cell(num_elems,length(t_nodes));
control_eval_nodes = cell(num_elems,length(t_nodes));

for i = 1:num_elems
  
    for q = 1:length(t_nodes)
        
        % storing as sparse matrices
        tmp = state_basis{i}(t_nodes(q));
        tmp(abs(tmp)<1e-9) = 0;         % numerical noise removal
        tmp = sparse(tmp);
        state_eval_nodes{i,q} = tmp;
        
        tmp = control_basis{i}(t_nodes(q));
        tmp(abs(tmp)<1e-9) = 0;         % numerical noise removal
        tmp = sparse(tmp);        
        control_eval_nodes{i,q} = tmp;
               
    end
    
end

end