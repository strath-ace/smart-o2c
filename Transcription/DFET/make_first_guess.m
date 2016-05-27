function x = make_first_guess(f,x_0,t_0,t_f,u,structure)

% check format of u (changes if using just fmincon or MACS)

if (size(u,1)*size(u,2))~=structure.num_controls*structure.num_elems*(structure.control_order+1)
   
    error('wrong number of control variables')
    
else % number of control variables is correct
    
    if size(u,2)>structure.num_controls % if ordering is not the expected matrix format but row vector
      
        % proper reshaping, good even for multi-elements, multi-controls
      
        u = reshape(u,structure.num_elems*(structure.control_order+1),structure.num_controls);
        for i = 1:structure.num_elems
            
            for j = 1:structure.num_controls
            
                %uu(1+(structure.control_order+1)*(i-1):(structure.control_order+1)*i,j) = u((1+((structure.control_order+1))*(i-1))*j:(structure.control_order+1)*i*j)';
            
            end
            
        end
        
    end
    
    %u = uu;
    
end

times = structure.in_nodes_state'*(t_f-t_0);
times = times(:);   % arranges times as a vector

% Clone Initial condition
x_temp = CloneIC(x_0,structure,times);

% Explicit Euler
%x_temp = ExpEulContr(f,x_0,u,structure,times);

% if max(max(abs(x_temp)))>1e15
%    
%     warning ('Explicit Euler found very large state values, initial guess might be VERY BAD. Try with more elements or change initial guess type')
%     
% end

% Implicit Euler
%x_temp = ImpEulContr(f,x_0,u,structure,times);

% Adaptive Explicit Euler
%x_temp = AdapExpEul(f,x_0,times);

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
    
    x(start_pos:end_pos) = reshape(u(1+(structure.control_order+1)*(i-1):(structure.control_order+1)*i,:),(structure.control_order+1)*structure.num_controls,1);
        
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