function [state_eval_nodes,control_eval_nodes,enodes] = eval_poly_on_nodes(num_elems,state_basis,state_nodes,control_basis,control_nodes)

% Evaluate polynomials of states and controls on the same given nodes

enodes = sort(unique([state_nodes,control_nodes]));

state_eval_nodes = cell(num_elems,length(enodes));
control_eval_nodes = cell(num_elems,length(enodes));

for i = 1:num_elems
  
    for q = 1:length(enodes)
        
        state_eval_nodes{i,q} = state_basis{i}(enodes(q));
        control_eval_nodes{i,q} = control_basis{i}(enodes(q));
               
    end
    
end

end