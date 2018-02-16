function x = make_first_guess(x_0,t_0,t_f,static,u,structure)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%
% Generates first guess. Different schemes implemented: cloning ic
% (simplest, less useful), DFET adding one element per time, and DFET all
% at once. The initial guess shouldn't automatically respect all
% constraints, obviously the more the better. DFET variations respect all
% constraints except for final equality constraints.
% NOTE: in the generation of the initial guess, static control parameters
% are for now considered fixed (and probably will stay so). So for instance
% inital/final time, initial conditions, etc should be supplied.
% DFET initializes states by cloning the IC. could be useful to randomly
% generate also the states, to have a more "free" exploration? Since it's
% basically little more than a straightforward integration, no convergence
% problems should arise...

% check if type is supported

type = structure.init_type;

if (~strcmp(type,'CloneIC'))&&(~strcmp(type,'DFET-all'))
    
    error('type must be CloneIC or DFET-all')
    
end

% check format of u (changes if using just fmincon or MACS)

if (size(u,1)*size(u,2))~=structure.num_controls*structure.num_elems*(structure.control_order+1)
    
    error('wrong number of control variables')
    
else % number of control variables is correct
    
    % proper reshaping, good even for multi-elements, multi-controls
    
    u = reshape(u,structure.num_controls*(structure.control_order+1),structure.num_elems);
    
end

if strcmp(type,'CloneIC')
    
    % Clone Initial condition
    
    times = structure.uniform_in_nodes_state'*(t_f-t_0);
    times = times(:);   % arranges times as a vector
    x_temp = CloneIC(x_0,structure,times);
    
    % rearranges solution as a vector (INCLUDING CONTROLS!!!!!)
    x = zeros(structure.num_elems*((structure.state_order+1)*structure.num_eqs+(structure.control_order+1)*structure.num_controls)+structure.num_free_states*(structure.DFET==1),1);
    
    start_pos = 1;
    
    for i = 1:structure.num_elems
        
        end_pos = start_pos+structure.num_eqs*(structure.state_order+1)-1;
        
        % all state coefficients for all equations of current element
        x(start_pos:end_pos) = reshape(x_temp(1+(structure.state_order+1)*(i-1):(structure.state_order+1)*i,:),(structure.state_order+1)*structure.num_eqs,1);
        
        % all control coefficients for all controls of current element
        start_pos = end_pos+1;
        end_pos = start_pos+structure.num_controls*(structure.control_order+1)-1;
        
        x(start_pos:end_pos) = u(:,i);%reshape(u(i,:),(structure.control_order+1)*structure.num_controls,1);
        
        start_pos = end_pos+1;
        
    end
    
    % inclusion of end condition for DFET, if it's free and thus an unknown
    
    if structure.DFET==1 && (structure.num_free_states>0)
        
        id_final = 1:structure.num_eqs;
        id_final = id_final(structure.free_final_states);
        
        x(start_pos:end) = x_temp(end,id_final);
        
    end
    
    % inclusion of static variables
    
    x_add = static;
    
    if ~structure.imposed_t0
        
        x_add = [x_add;t_0];
        
    end
    
    x_add = [x_add; x_0(~structure.imposed_initial_states)];
    x = [x_add;x];
    
    if ~structure.imposed_tf
        
        x = [x; t_f];
        
    end
    
else
    
    x = DFETIC2(x_0,t_0,t_f,static,u,structure);
    
end

end