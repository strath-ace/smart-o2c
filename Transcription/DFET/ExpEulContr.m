function y = ExpEulContr(f,x_0,u_state,structure,t)

% Integrates with Explicit Euler over PRESCRIBED vector of times t

num_times = length(t)-1*(structure.DFET==1)*(structure.state_order==0);
num_free_end_states = sum(structure.free_final_states);
coeffs = structure.state_order+1;

if size(x_0,1)>1
   
    x_0 = x_0';
    
end

y = zeros(num_times+(num_free_end_states>0)*(structure.DFET==1),structure.num_eqs);
y(1,:) = x_0;

for i =2:num_times

    % locate in which element is located current timestep
    curr_elem = ceil(i/(coeffs+1*(structure.state_order==0)));
       
    % extrapolate correct subvector
    this_u = u_state(1+(structure.control_order+1)*structure.num_controls*(curr_elem-1):(structure.control_order+1)*curr_elem*structure.num_controls,:);
    
    % scale time between to element-time [-1 1]
    tmax = structure.in_nodes_state(curr_elem,end)*t(end);
    tmin = structure.in_nodes_state(curr_elem,1)*t(end);
    
    norm_time = 2*(t(i)-(tmax+tmin)/2)/(tmax-tmin);
        
    % evaluate u
    u_now = structure.control_basis{curr_elem}(norm_time)*this_u;
    
    % update
    y(i,:) = y(i-1,:)+f(y(i-1,:),u_now,t(i))'*(t(i)-t(i-1));
    
end

%this has to be improved/specialized
if structure.DFET==1 && (num_free_end_states>0) 

    y(end,:) = y(end-1,:);
    
end

end