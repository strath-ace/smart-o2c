function x = update_guess(u,x_temp,structure)

x = zeros(structure.num_elems*((structure.state_order+1)*structure.num_eqs+(structure.control_order+1)*structure.num_controls)+structure.num_free_states*(structure.DFET==1),1);

start_pos = 1;

for i = 1:structure.num_elems
    
    end_pos = start_pos+structure.num_eqs*(structure.state_order+1)-1;
    
    % all state coefficients for all equations of current element
    x(start_pos:end_pos) = reshape(x_temp(:,i),(structure.state_order+1)*structure.num_eqs,1);
    
    % all control coefficients for all controls of current element
    start_pos = end_pos+1;
    end_pos = start_pos+structure.num_controls*(structure.control_order+1)-1;
    
    x(start_pos:end_pos) = u(1+(structure.control_order+1)*structure.num_controls*(i-1):(structure.control_order+1)*structure.num_controls*i);
    
    start_pos = end_pos+1;
    
end

% inclusion of end condition for DFET, if it's free and thus an unknown

if structure.DFET==1 && (structure.num_free_states>0)
    
    id_final = 1:structure.num_eqs;
    id_final = id_final(structure.free_final_states);
    
    x(start_pos:end) = x_temp(end,id_final);
       
end

%x = x(:);

end