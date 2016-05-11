function y = ImpEulContr(f,x_0,u_state,structure,t)

% Integrates with Implicit Euler over PRESCRIBED vector of times t

num_times = length(t)-1*(structure.DFET==1)*(structure.state_order==0);
num_free_end_states = sum(structure.free_final_states);
coeffs = structure.state_order+1;

if size(x_0,1)>1
   
    x_0 = x_0';
    
end

options = optimset('Display','iter');

y = zeros(num_times+(num_free_end_states>0)*(structure.DFET==1),structure.num_eqs);
y(1,:) = x_0;

for i =2:num_times

    % locate in which element is located current timestep
    curr_elem = ceil(i/(coeffs+1*(structure.state_order==0)));
       
    % extrapolate correct subvector
    this_u = u_state(1+(structure.control_order+1)*structure.num_controls*(curr_elem-1):(structure.control_order+1)*curr_elem*structure.num_controls);
    
    % scale time between to element-time [-1 1]
    tmax = structure.in_nodes_state(curr_elem,end);
    tmin = structure.in_nodes_state(curr_elem,1);
    
    norm_time = 2*(t(i)-(tmax+tmin)/2)/(tmax-tmin);
        
    % evaluate u
    u_now = structure.control_basis{curr_elem}(norm_time)*this_u;
    
    % update
    ypred =  y(i-1,:);%+f(y(i-1,:),u_now,t(i))'*(t(i)-t(i-1));
    sol = fsolve(@(x) funToSolve(ypred',u_now',t(i),y(i-1,:)',f,t(i)-t(i-1)),ypred',options);

    y(i,:) = sol';
    %y(i,:) = y(i-1,:)+f(y(i-1,:),u_now,t(i))'*(t(i)-t(i-1));
    
end

%this has to be improved/specialized
if structure.DFET==1 && (num_free_end_states>0) 

    y(end,:) = y(end-1,:);
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The function to solve 
function residual = funToSolve(x,u,t,xo,MyFunc,dt)
  residual=xo+feval(MyFunc,x,u,t)*dt-x;
return

end