function [lbv,ubv,state_vars,control_vars] = transcribe_bounds(state_bounds,control_bounds,structure)

%% Input check

if size(state_bounds,1)~=structure.num_eqs
   
    error('state_bounds must have as many rows as equations')
    
end

if size(state_bounds,2)~=2
   
    error('state_bounds must have 2 columns per row [lb ub]')
    
end

if size(control_bounds,1)~=structure.num_controls
   
    error('control_bounds must have as many rows as controls')
    
end

if size(control_bounds,2)~=2
   
    error('control_bounds must have 2 columns per row [lb ub]')
    
end

%% Creation of ordered lower/upper bound vector

pad = zeros(structure.num_eqs*(structure.state_order+1)+structure.num_controls*(structure.control_order+1),2);
state_pad = pad(:,1);
control_pad = pad(:,1);

% create padding matrix, will simply be cloned for all elements

for i = 1:structure.num_eqs
    
    pad(1+(structure.state_order+1)*(i-1):(structure.state_order+1)*(i),:) = repmat(state_bounds(i,:),structure.state_order+1,1);
    state_pad(1+(structure.state_order+1)*(i-1):(structure.state_order+1)*(i),:) = 1;
    
end

for i = 1:structure.num_controls

    pad((structure.state_order+1)*structure.num_eqs+1+(structure.control_order+1)*(i-1):(structure.state_order+1)*structure.num_eqs+(structure.control_order+1)*(i),:) = repmat(control_bounds(i,:),structure.control_order+1,1);
    control_pad((structure.state_order+1)*structure.num_eqs+1+(structure.control_order+1)*(i-1):(structure.state_order+1)*structure.num_eqs+(structure.control_order+1)*(i),:) = 1;
    
end

lbv = repmat(pad(:,1),structure.num_elems,1);
ubv = repmat(pad(:,2),structure.num_elems,1);
state_vars = repmat(state_pad(:,1),structure.num_elems,1);
control_vars = repmat(control_pad(:,1),structure.num_elems,1);

if structure.DFET==1
   
    for i = 1:length(structure.imposed_final_states)
    
        if ~(structure.imposed_final_states(i))
        
            lbv = [lbv; pad((structure.state_order+1)*i,1)];
            ubv = [ubv; pad((structure.state_order+1)*i,2)];
            state_vars = [state_vars;state_pad((structure.state_order+1),1)];
            control_vars = [control_vars;0*state_pad((structure.state_order+1),1)];
            
        end
        
    end
    
end

end