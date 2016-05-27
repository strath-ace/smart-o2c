function y = CloneIC(x_0,structure,t)

% Integrates with Implicit Euler over PRESCRIBED vector of times t

num_times = length(t)-1*(structure.DFET==1)*(structure.state_order==0);
num_free_end_states = sum(structure.free_final_states);

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