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

mask = zeros((structure.control_order+1)*size(u_state,2),size(u_state,2));

for j = 1:size(u_state,2)

    mask(1+(structure.control_order+1)*(j-1):(structure.control_order+1)*j,j) = ones(structure.control_order+1,1);
    
end

for i =2:num_times

    % locate in which element is located current timestep
    curr_elem = ceil(i/(coeffs+1*(structure.state_order==0)));
       
    % extrapolate correct subvector
    this_u = u_state(1+(structure.control_order+1)*(curr_elem-1):(structure.control_order+1)*curr_elem,:);
    
    % correct size of this_u when more than one control parameter is
    % present
    
    this_u = mask.*repmat(this_u,structure.num_controls,1);
    
    % scale time between to element-time [-1 1]
    tmax = structure.in_nodes_state(curr_elem,end)*t(end);
    tmin = structure.in_nodes_state(curr_elem,1)*t(end);
    
    % old version, t(i), probably wrong
    %norm_time = 2*(t(i)-(tmax+tmin)/2)/(tmax-tmin);
        
    norm_time = 2*(t(i-1)-(tmax+tmin)/2)/(tmax-tmin);
    
    % evaluate u
    u_now = structure.control_basis{curr_elem}(norm_time)*this_u;
    
    % update
    % old version, t(i), probably wrong
    %y(i,:) = y(i-1,:)+f(y(i-1,:),u_now,t(i))'*(t(i)-t(i-1));
   
    y(i,:) = y(i-1,:)+f(y(i-1,:),u_now,t(i-1))'*(t(i)-t(i-1));
    
end

%this has to be improved/specialized
if structure.DFET==1 && (num_free_end_states>0) 

    y(end,:) = y(end-1,:);
    
end

end