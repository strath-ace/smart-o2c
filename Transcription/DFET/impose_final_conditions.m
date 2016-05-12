function structure = impose_final_conditions(structure, imposed_final_states)

% Inputs should already be checked at this point

% Continuous Galerkin adds matching equations, Discontinuous Galerkin adds
% unconstrained final states as unknowns!!!!!!!

num_imposed_states = sum (imposed_final_states);

id_imposed_states = 1:length(imposed_final_states);
id_imposed_states = id_imposed_states(imposed_final_states==1);    % this vector basically contains the id of the "equations" where the final states are imposed

num_free_states = sum (~imposed_final_states);

id_free_states = 1:length(imposed_final_states);
id_free_states = id_free_states(imposed_final_states==0);    % this vector basically contains the id of the "equations" where the final states are imposed


if (length(id_imposed_states)~= num_imposed_states)
    
    error('Something is wrong...');
    
end

if structure.DFET == 0
    
    % initialise augmented matrix
    M_aug = zeros(size(structure.M,1)+num_imposed_states,size(structure.M,2));
    
    % copy all but last element
    M_aug(1:(structure.state_order+1)*structure.num_eqs*(structure.num_elems-1),1:(structure.state_order+1)*structure.num_eqs*(structure.num_elems-1)) = structure.M(1:(structure.state_order+1)*structure.num_eqs*(structure.num_elems-1),1:(structure.state_order+1)*structure.num_eqs*(structure.num_elems-1));
    
    row_condition = structure.state_basis{end}(1);
    
    if structure.num_elems == 1 % has to be done slightly differently since there's no "backwards match condition" to impose...
        
        xstart_position_aug = 1;
        ystart_position_aug = 1;
        
        xstart_position = 1;
        ystart_position = 1;
        
        for i = 1:structure.num_eqs
            
            xend_position = xstart_position + (structure.state_order+1)-1;
            yend_position = xend_position;
            
            M_temp = structure.M(ystart_position:yend_position,xstart_position:xend_position);
            
            if imposed_final_states(i)==1
                
                this_row_condition = row_condition(i,(i-1)*(structure.state_order+1)+1:i*(structure.state_order+1));
                
                M_temp = [M_temp; this_row_condition];
                
            end
            
            xend_position_aug = xstart_position_aug + size(M_temp,2)-1;
            yend_position_aug = ystart_position_aug + size(M_temp,1)-1;
            
            
            % Update augmented matrix
            
            M_aug(ystart_position_aug:yend_position_aug,xstart_position_aug:xend_position_aug) = M_temp;
            
            % initialise start_position for next loop
            
            xstart_position_aug = xstart_position_aug+structure.state_order+1;
            ystart_position_aug = ystart_position_aug+structure.state_order+1+1*(imposed_final_states(i)==1);
            
            xstart_position = xstart_position+(structure.state_order+1);
            ystart_position = ystart_position+(structure.state_order+1);
            
        end
        
    else
        
        ystart_position_aug = (structure.state_order+1)*(structure.num_elems-1)*structure.num_eqs+1;
        xstart_position_aug = ystart_position_aug-(structure.state_order+1)*(structure.num_eqs);
        
        xstart_position = xstart_position_aug;
        ystart_position = ystart_position_aug;
        
        
        for i = 1:structure.num_eqs
            
            % identify position within the matrices
            
            xend_position = xstart_position + (structure.state_order+1)*(structure.num_eqs+1)-1;
            yend_position = ystart_position + (structure.state_order);
            
            % submatrix to be copied, contains submatrices of all previous
            % equations with NO imposed end condition, and current equation which
            % has imposed condition plus the imposed condition
            
            M_temp = structure.M(ystart_position:yend_position,xstart_position:xend_position);
            
            if (imposed_final_states(i)==1)
                
                % condition to be imposed
                this_row_condition = [zeros(1,((structure.state_order+1)*(structure.num_eqs))) row_condition(i,(i-1)*(structure.state_order+1)+1:i*(structure.state_order+1))];
                M_temp = [M_temp; this_row_condition];
                
            end
            
            xend_position_aug = xstart_position_aug + size(M_temp,2)-1;
            yend_position_aug = ystart_position_aug + size(M_temp,1)-1;
            
            % Update augmented matrix
            
            M_aug(ystart_position_aug:yend_position_aug,xstart_position_aug:xend_position_aug) = M_temp;
            
            % initialise start_position for next loop
            
            xstart_position_aug = xstart_position_aug+structure.state_order+1;
            ystart_position_aug = ystart_position_aug+structure.state_order+1+1*(imposed_final_states(i)==1);
            
            xstart_position = xstart_position+(structure.state_order+1);
            ystart_position = ystart_position+(structure.state_order+1);
            
        end
        
    end
    
else
    
    % initialize augmented Matrix
    M_aug = zeros(size(structure.M,1),num_free_states+size(structure.M,2));
    
    % copy all existing matrix in the augmented one
    
    M_aug(:,1:end-num_free_states) = structure.M;
    
    % projection of test function on right edge of element
    vals = structure.test_basis{end}(1)';
    
    if structure.num_elems==1
        
        for j = 1:num_free_states
            
            % include free end conditions as additional variables
            
            M_aug(end-size(vals,1)+1:end,size(structure.M,2)+j) = -vals(:,id_free_states(j));% -structure.tes;
            
        end
        
    else
        
        starty = size(M_aug,1)-structure.num_eqs*structure.test_order;
        startx = size(structure.M,2);
        
        if num_free_states>0
            
            for j = 1:length(imposed_final_states)
                
                startx = startx + 1*(imposed_final_states(j)==0);
                M_aug(starty+j*structure.test_order,startx) = -1*(imposed_final_states(j)==0);
                
            end
            
        end
        
    end
    
end

%% Cleanup of should-be-zero elements (truncating entries < 1e-15 to 0)

M_aug(abs(M_aug)<1e-10)=0;

structure.M = M_aug;
structure.imposed_final_states = imposed_final_states;
structure.num_imposed_states = sum(imposed_final_states);
structure.free_final_states = ~imposed_final_states;
structure.num_free_states = sum(structure.free_final_states);


end