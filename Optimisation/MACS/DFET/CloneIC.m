function y = CloneIC(f,x_0,u_state,structure,t)

% Integrates with Implicit Euler over PRESCRIBED vector of times t

num_times = length(t)-1*(structure.DFET==1)*(structure.state_order==0);
num_free_end_states = sum(structure.free_final_states);
coeffs = structure.state_order+1;

if size(x_0,1)>1
   
    x_0 = x_0';
    
end

y = zeros(num_times+(num_free_end_states>0)*(structure.DFET==1),structure.num_eqs);
y(1,:) = x_0;

for i =2:num_times

    y(i,:) = x_0;
    
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