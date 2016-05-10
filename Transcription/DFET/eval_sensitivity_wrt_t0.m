function F = eval_sensitivity_wrt_t0(f,structure,t_0,t_f,loc_state,loc_u,dft)

% evals sensitivity of x(t_f) (intended as both the states and the
% controls) wrt t_0. BEWARE, IT'S USEFUL ONLY FOR THE GRADIENT OF THE
% OBJECTIVE FUNCTIONS EVALUATED AT THE FINAL TIME (g(t_f))

% computation of integrals. These DO influence the computation of
% derivatives. In particular, these integrals change when t_0, t_f and the
% states change. x_0 and x_f DO NOT change the integrals, instead.

F = zeros(size(structure.M,1),1);

ee1 = structure.els(:,1)*(t_f-t_0);
de1 = -structure.els(:,1);
ee2 = structure.els(:,2)*(t_f-t_0);
de2 = -structure.els(:,2);

tmpp = cell(structure.num_elems,1);

for i = 1:structure.num_elems
    
    tmpp{i} = 0;
    
    e1 = ee1(i);
    d1 = de1(i);
    e2 = ee2(i);
    d2 = de2(i);
    
    % maybe it's better to avoid reshape, saves a for loop afterwards for
    % the reconstruction of complete f
    
    for j = 1:length(structure.integr_nodes)
        
        tmpp{i} = tmpp{i} + reshape(structure.integr_weights(j)*structure.test_eval{i,j}'*f(structure.state_eval{i,j}*loc_state(:,i),structure.control_eval{i,j}*loc_u(:,i),(e2+e1)/2+(e2-e1)/2*structure.integr_nodes(j)).*(d2-d1)/2,structure.test_order+1,structure.num_eqs);
        tmpp{i} = tmpp{i} + reshape(structure.integr_weights(j)*structure.test_eval{i,j}'*dft(structure.state_eval{i,j}*loc_state(:,i),structure.control_eval{i,j}*loc_u(:,i),(e2+e1)/2+(e2-e1)/2*structure.integr_nodes(j)).*((d2+d1)/2+(d2-d1)/2*structure.integr_nodes(j)).*(e2-e1)/2,structure.test_order+1,structure.num_eqs);
        
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
            F_temp(1) = 0;
            
            if (structure.imposed_final_states(q)==1)
                
                F_temp = [F_temp;0];               % matching condition <- FOR THE FET THIS IS THE ONLY POINT WHERE x_f ENTERS, AND IS IN A LINEAR FASHION!!
                
            end
            
            yend_aug = ystart_aug + length(F_temp)-1;
            
            F(ystart_aug:yend_aug) = F_temp;
            
            ystart_aug = yend_aug+1;
            
        end
        
    else
        
        ystart_aug = 1;        
        
        for q = 1:structure.num_eqs
            
            F_temp = -tmpp{end}(:,q);
                        
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
            
            F(1+(structure.state_order+1)*(q-1)) = 0;
            
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
                
                F_temp = [F_temp;0];
                
            end
            
            yend_aug = ystart_aug + length(F_temp)-1;
            
            F(ystart_aug:yend_aug) = F_temp;
            
            ystart_aug = yend_aug+1;
            
        end
        
    else
        
        % first element, contains initial conditions
        
        ystart_aug = 1;
        
        %evaluation of test basis on -1, to introduce initial conditions
        
        for q = 1:structure.num_eqs
            
            F_temp = -tmpp{1}(:,q);
            
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
                
        start2 = start2+structure.test_order+1*(structure.num_elems==2);
        
        for q = 1:structure.num_eqs
            
            F_temp = -tmpp{end}(:,q);
                        
            yend_aug = ystart_aug + length(F_temp)-2;
            
            F(ystart_aug:yend_aug) = F_temp(2:end);
            F(start2) = F(start2)+F_temp(1);
            
            ystart_aug = yend_aug+1;
            start2 = start2+structure.test_order;
            
        end
        
    end
    
end

end