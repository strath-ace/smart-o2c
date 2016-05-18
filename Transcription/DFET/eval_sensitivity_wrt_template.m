function ddx0 = eval_sensitivity_wrt_template(f,structure,x,x_0,x_f,ee1,ee2,loc_state,loc_u)

% evals sensitivity of x wrt TEMPLATE

% computation of integrals. These DO influence the computation of
% derivatives. In particular, these integrals change when t_0, t_f and the
% states change. x_0 and x_f DO NOT change the integrals, instead.

for i = 1:structure.num_elems
    
    e1 = ee1(i);
    e2 = ee2(i);
    
    % maybe it's better to avoid reshape, saves a for loop afterwards for
    % the reconstruction of complete f
    
    for j = 1:length(structure.integr_nodes)
        
        tmpp{i} = tmpp{i} + reshape(structure.integr_weights(j)*structure.test_eval{i,j}'*f(structure.state_eval{i,j}*loc_state(:,i),structure.control_eval{i,j}*loc_u(:,i),(e2+e1)/2+(e2-e1)/2*structure.integr_nodes(j)).*(e2-e1)/2,structure.test_order+1,structure.num_eqs);
        
    end
    
end

% reordering of the solution, considering initial, matching and final
% conditions. Here all the variations are involved (even wrt x_0 and x_f).
% It's better to leave the order of the elements unchanged and extract the
% derivatives wrt x and u afterwards

if structure.num_elems == 1
    
    if structure.DFET==0
        
        ystart_aug = 1;
        
        for q = 1:structure.num_eqs
            
            F_temp = tmpp{end}(:,q);
            F_temp(1) = x_0(q);      % matching condition <- FOR THE FET THIS IS THE ONLY POINT WHERE x_0 ENTERS, AND IS IN A LINEAR FASHION!!
            
            if (structure.imposed_final_states(q)==1)
                
                F_temp = [F_temp;x_f(q)];               % matching condition <- FOR THE FET THIS IS THE ONLY POINT WHERE x_f ENTERS, AND IS IN A LINEAR FASHION!!
                
            end
            
            yend_aug = ystart_aug + length(F_temp)-1;
            
            F(ystart_aug:yend_aug) = F_temp;
            
            ystart_aug = yend_aug+1;
            
        end
        
    else
        
        ystart_aug = 1;
        
        %evaluation of test basis on -1, to introduce initial conditions
        valsl = structure.test_basis{end}(-1)';
        %evaluation of test basis on 1, to introduce known final conditions
        valsr = structure.test_basis{end}(1)';
        
        for q = 1:structure.num_eqs
            
            F_temp = -tmpp{end}(:,q)-valsl(1+(structure.test_order+1)*(q-1):(structure.test_order+1)*q,q)*x_0(q);
            
            if (structure.imposed_final_states(q)==1)
                
                F_temp = F_temp+valsr(1+(structure.test_order+1)*(q-1):(structure.test_order+1)*q,q)*x_f(q);
                
            end
            
            yend_aug = ystart_aug + length(F_temp)-1;
            
            F(ystart_aug:yend_aug) = F_temp;
            
            ystart_aug = yend_aug+1;
            
        end
        
    end
    
else
    
    if structure.DFET==0 % sarebbe piu' sano diversificare le eval_constraints, tanto DFET non puo' mai cambiare ma questo if viene valutato ogni volta...
        
        % copy solutions, up to the last element, where modifications of
        % ordering occour
        
        for i=1:structure.num_elems-1
            
            F(1+(structure.test_order+1)*structure.num_eqs*(i-1):(structure.test_order+1)*structure.num_eqs*(i)) = reshape(tmpp{i},(structure.test_order+1)*structure.num_eqs,1);
            
        end
        
        
        % first element, impose initial conditions
        
        for q = 1:structure.num_eqs
            
            F(1+(structure.state_order+1)*(q-1)) = x_0(q);
            
        end
        
        % intermediate elements, impose matching
        
        for i = 2:structure.num_elems-1
            
            for q = 1:structure.num_eqs
                
                F(1+(q-1)*(structure.state_order+1)+(structure.state_order+1)*structure.num_eqs*(i-1)) = 0;
                
            end
            
        end
        
        % last element, impose matching, reorder equations and impose final
        % state
        
        ystart_aug = (structure.state_order+1)*structure.num_eqs*(structure.num_elems-1)+1;
        
        for q = 1:structure.num_eqs
            
            F_temp = tmpp{end}(:,q);
            F_temp(1) = 0;      % matching condition
            
            if (structure.imposed_final_states(q)==1)
                
                F_temp = [F_temp;x_f(q)];
                
            end
            
            yend_aug = ystart_aug + length(F_temp)-1;
            
            F(ystart_aug:yend_aug) = F_temp;
            
            ystart_aug = yend_aug+1;
            
        end
        
    else
        
        % first element, contains initial conditions
        
        ystart_aug = 1;
        
        %evaluation of test basis on -1, to introduce initial conditions
        valsl = structure.test_basis{1}(-1)';
        
        for q = 1:structure.num_eqs
            
            F_temp = -tmpp{1}(:,q)-valsl(1+(structure.test_order+1)*(q-1):(structure.test_order+1)*q,q)*x_0(q);
            
            yend_aug = ystart_aug + length(F_temp)-1;
            
            F(ystart_aug:yend_aug) = F_temp;
            
            ystart_aug = yend_aug+1;
            
        end
        
        % central elements, copy all integrals except for the first one,
        % that has to be added to the corresponding of the previous element
        % (messy)
        
        start2 = 0;
        
        if structure.num_elems>2
            
            for i = 2:structure.num_elems-1
                
                for q = 1:structure.num_eqs
                    
                    start2 = start2+structure.test_order+1*(i==2);
                    
                    F_temp = -tmpp{i}(:,q);
                    
                    yend_aug = ystart_aug + length(F_temp)-2;
                    
                    F(ystart_aug:yend_aug) = F_temp(2:end);
                    F(start2) = F(start2)+F_temp(1);
                    
                    ystart_aug = yend_aug+1;
                    
                end
                
            end
            
        end
        
        % final element
        
        valsr = structure.test_basis{end}(1)';
        
        start2 = start2+structure.test_order+1*(structure.num_elems==2);
        
        for q = 1:structure.num_eqs
            
            F_temp = -tmpp{end}(:,q);
            
            if (structure.imposed_final_states(q)==1)
                
                F_temp = F_temp+valsr(1+(structure.test_order+1)*(q-1):(structure.test_order+1)*q,q)*x_f(q);
                
            end
            
            yend_aug = ystart_aug + length(F_temp)-2;
            
            F(ystart_aug:yend_aug) = F_temp(2:end);
            F(start2) = F(start2)+F_temp(1);
            
            ystart_aug = yend_aug+1;
            start2 = start2+structure.test_order;
            
        end
        
    end
    
end

x_only = loc_state(:);

if structure.DFET==1
    
    x_only = [x_only;x(end-structure.num_free_states+1:end)'];
    
end

F = structure.M*x_only-F;

end