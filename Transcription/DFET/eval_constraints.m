function [F,J] = eval_constraints(f,structure,x,x_0,x_f,t_0,t_f,els,calc_jac,varargin)

dfx = varargin{1};
dfu = varargin{2};

if xor(isempty(dfx),isempty(dfu))
    
    error('Either give analytic Jacobian of f with respect to x and u, or give none')
    
end

%% Computation of the RHS

%F = zeros(structure.num_elems*structure.num_eqs*(structure.test_order+1)+sum(structure.imposed_final_states)+structure.num_free_states*(structure.DFET==1),1);

F = zeros(size(structure.M,1),1);

ee1 = els(:,1)*(t_f-t_0);%structure.els(:,1)*(t_f-t_0);
ee2 = els(:,2)*(t_f-t_0);%structure.els(:,2)*(t_f-t_0);

loc_state = zeros((structure.state_order+1)*structure.num_eqs,structure.num_elems);
loc_u = zeros((structure.control_order+1)*structure.num_controls,structure.num_elems);

tmpp = cell(structure.num_elems,1);

startx = 1;

for i=1:structure.num_elems
    
    tmpp{i} = 0;%zeros(structure.state_order+1,structure.num_eqs);
    
    endx = startx+(structure.state_order+1)*structure.num_eqs-1;
    loc_state(:,i) = x(startx:endx);
    
    startx = endx+1*(structure.num_controls>0);
    endx = startx+(structure.control_order+1)*structure.num_controls-1*(structure.num_controls>0);
    
    loc_u(:,i) = x(startx:endx);
    
    startx = endx+1;
end

% compute integrals elementwise (can be made PARALLEL, although benefits depend widely upon the actual amount of work done by each node)

for i = 1:structure.num_elems
    
    e1 = ee1(i);
    e2 = ee2(i);
    
    % maybe it's better to avoid reshape, saves a for loop afterwards for
    % the reconstruction of complete f
    
    for j = 1:length(structure.integr_nodes)
        
        tmpp{i} = tmpp{i} + reshape(structure.integr_weights(j)*structure.test_eval{i,j}'*f(structure.state_eval{i,j}*loc_state(:,i),structure.control_eval{i,j}*loc_u(:,i),(e2+e1)/2+(e2-e1)/2*structure.integr_nodes(j)).*(e2-e1)/2,structure.test_order+1,structure.num_eqs);
        
    end
    
end

% reordering of the solution, considering initial, matching and final conditions

if structure.num_elems == 1
    
    if structure.DFET==0
        
        ystart_aug = 1;
        
        for q = 1:structure.num_eqs
            
            F_temp = tmpp{end}(:,q);
            F_temp(1) = x_0(q);      % matching condition
            
            if (structure.imposed_final_states(q)==1)
                
                F_temp = [F_temp;x_f(q)];
                
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
    
    if size(x,1)==1
       
        xx = x(end-structure.num_free_states+1:end)';
        
    else
        
        xx = x(end-structure.num_free_states+1:end);
        
    end
    
    x_only = [x_only;xx];
    
end

F = structure.M*x_only-F;

%% Computation of the Jacobian in the projected state space

J = zeros(length(F),length(x));

if calc_jac
    
    if isempty(dfx)
        
        %% FINITE DIFFERENCES
                
        % first order "forward" finite difference approach with fixed step
        % of 1e-6. Would like to use backwards differences when variable is
        % at upper bound, but would require to know the bounds here. Second
        % order accurate finite differences could also be appealing,
        % requiring to average forward and backward differences and
        % eventually using ad adapted stencil for the boundaries, but
        % probably the additional computational cost does not pay off.        
        for i = 1:length(x)
            
            x_v = x;
            x_v(i) = x_v(i)+0.000001;
            F_temp = eval_constraints(f,structure,x_v,x_0,x_f,t_0,t_f,els,0,dfx,dfu);
            J(:,i) = (F_temp-F)/0.000001;
            
        end
        
    else
        
        %J = zeros(length(F),structure.num_elems*(structure.num_eqs*(structure.state_order+1)+structure.num_controls*(structure.control_order+1)));
        
        %% ANALYTIC JACOBIAN
        
        % compute Jacobians wrt state variables in projected space, elementwise
        
        tmpJx = cell(structure.num_elems,1);
        
        for i = 1:structure.num_elems
            
            tmpJx{i} = 0;
            
            for j = 1:length(structure.integr_nodes)
                
                tmpJx{i}= tmpJx{i} + structure.integr_weights(j)*structure.test_eval{i,j}'*dfx(structure.state_eval{i,j}*loc_state(:,i),structure.control_eval{i,j}*loc_u(:,i),(e2+e1)/2+(e2-e1)/2*structure.integr_nodes(j))*structure.state_eval{i,j}.*(e2-e1)/2;
                
            end
            
            if structure.DFET ==0
                % remove first row for each equation of each element (matching conditions will be applied)
                for q = 1:structure.num_eqs
                    
                    tmpJx{i}(1+(structure.state_order+1)*(q-1),:) = 0*tmpJx{i}(1+(structure.state_order+1)*(q-1),:);
                    
                end
                
            end
            
        end
        
        % compute Jacobians wrt control variables (if any!!) in projected space, elementwise
        
        if structure.num_controls>0
            
            tmpJu = cell(structure.num_elems,1);
            
            for i = 1:structure.num_elems
                
                tmpJu{i} = 0;
                
                for j = 1:length(structure.integr_nodes)
                    
                    tmpJu{i}= tmpJu{i} + structure.integr_weights(j)*structure.test_eval{i,j}'*dfu(structure.state_eval{i,j}*loc_state(:,i),structure.control_eval{i,j}*loc_u(:,i),(e2+e1)/2+(e2-e1)/2*structure.integr_nodes(j))*structure.control_eval{i,j}.*(e2-e1)/2;
                    
                end
                
                if structure.DFET ==0
                    % remove first row for each element (matching conditions will be applied)
                    for q = 1:structure.num_eqs
                        
                        tmpJu{i}(1+(structure.state_order+1)*(q-1),:) = 0*tmpJu{i}(1+(structure.state_order+1)*(q-1),:);
                        
                    end
                    
                end
                
            end
            
        end
        
        % reordering of the solution, considering initial, matching and final conditions
        
        if structure.num_elems==1
            
            if structure.DFET==0
                
                ystart_aug = 1;
                
                for q = 1:structure.num_eqs
                    
                    Jfxq = structure.M(ystart_aug:ystart_aug+(structure.state_order),:)-tmpJx{end}(1+(structure.state_order+1)*(q-1):(structure.state_order+1)*q,:);
                    %Jfxq(1,:) = zeros(1,size(Jfxq,2));      % matching (initial) condition
                    %Jfxq(1,1+(structure.state_order+1)*(q-1)) = 1;
                    
                    if (structure.imposed_final_states(q)==1)
                        
                        endval = structure.state_basis{end}(1);
                        endval = endval(q,:);
                        Jfxq = [Jfxq; endval];
                        %Jfxq(end,(structure.state_order+1)*q) = 1;
                        
                    end
                    
                    % include Jacobian of controls, if present
                    
                    if structure.num_controls>0
                        
                        Jfuq = -tmpJu{end}(1+(structure.state_order+1)*(q-1):(structure.state_order+1)*q,:);
                        %Jfuq(1,:) = zeros(1,size(Jfuq,2));
                        
                        if (structure.imposed_final_states(q)==1)
                            
                            Jfuq = [Jfuq; zeros(1,size(Jfuq,2))];
                            
                        end
                        
                        
                    end
                    
                    if structure.num_controls>0
                        
                        J_temp = [Jfxq Jfuq];
                        
                    else
                        
                        J_temp = Jfxq;
                        
                    end
                    
                    
                    yend_aug = ystart_aug + size(J_temp,1)-1;
                    
                    J(ystart_aug:yend_aug,:) = J_temp;
                    
                    ystart_aug = yend_aug+1;
                    
                    
                end
                
            else
                
                ystart_aug = 1;
                
                for q = 1:structure.num_eqs
                    
                    Jfxq = structure.M(ystart_aug:ystart_aug+(structure.test_order),:)+[tmpJx{1}(1+(structure.test_order+1)*(q-1):(structure.test_order+1)*q,:) zeros(structure.test_order+1,structure.num_free_states)];
                    
                    % include Jacobian of controls, if present
                    
                    if structure.num_controls>0
                        
                        Jfuq = tmpJu{end}(1+(structure.test_order+1)*(q-1):(structure.test_order+1)*q,:);
                        
                    end
                    
                    if structure.num_controls>0
                        
                        if structure.num_free_states==0
                            
                            J_temp = [Jfxq Jfuq];
                            
                        else
                            
                            J_temp = [Jfxq(:,1:end-structure.num_free_states) Jfuq Jfxq(:,end-structure.num_free_states+1:end)];
                            
                        end
                        
                    else
                        
                        J_temp = Jfxq;
                        
                    end
                    
                    yend_aug = ystart_aug + size(J_temp,1)-1;
                    
                    J(ystart_aug:yend_aug,:) = J_temp;
                    
                    ystart_aug = yend_aug+1;
                    
                end
                
            end
            
        else
            
            if structure.DFET==0
                
                startx = 1;
                starty = 1;
                
                % copy solution for all but last elements
                for i = 1:structure.num_elems-1
                    
                    Jtmp = structure.M(1+(structure.state_order+1)*(i-1)*structure.num_eqs:(structure.state_order+1)*i*structure.num_eqs,1+(structure.state_order+1)*(i-1)*structure.num_eqs:(structure.state_order+1)*i*structure.num_eqs) - tmpJx{i};
                    
                    % add Jacobian of controls, if present
                    if structure.num_controls>0
                        
                        Jtmp = [Jtmp tmpJu{i}];
                        
                    end
                    
                    endx = startx+size(Jtmp,2)-1;
                    endy = starty+size(Jtmp,1)-1;
                    
                    J(starty:endy,startx:endx) = Jtmp;
                    
                    starty = endy+1;
                    startx = endx+1;
                    
                end
                
                % apply matching condition to elements 2 to end-1
                
                startx = 1;
                starty = 1+(structure.state_order+1)*structure.num_eqs;
                
                for i = 2:structure.num_elems-1
                    
                    endval = -structure.state_basis{i-1}(1);
                    
                    for q = 1:structure.num_eqs
                        
                        thisval = endval(q,:);
                        
                        J(starty,startx:startx+length(thisval)-1) = thisval;
                        
                        %startx = startx+(structure.state_order+1);
                        starty = starty+(structure.state_order+1);
                        
                    end
                    
                    startx = startx+(structure.state_order+1)*structure.num_eqs+(structure.control_order+1)*structure.num_controls;
                    
                end
                
                % apply matching and final conditions to last element
                
                startx_aug = endx+1;
                starty_aug = endy+1;
                
                starty = (structure.num_elems-1)*structure.num_eqs*(structure.state_order+1)+1;
                startx = starty;
                
                endval = structure.state_basis{end}(1);
                prevval = -structure.state_basis{end-1}(1);
                
                for q = 1:structure.num_eqs
                    
                    Jfxq = structure.M(starty:starty+(structure.state_order),startx:end)-tmpJx{end}(1+(structure.state_order+1)*(q-1):(structure.state_order+1)*q,:);
                    %Jfxq(1,:) = zeros(1,size(Jfxq,2));      % matching (initial) condition
                    %Jfxq(1,1+(structure.state_order+1)*(q-1)) = 1;
                    
                    if (structure.imposed_final_states(q)==1)
                        
                        thisval = endval(q,:);
                        Jfxq = [Jfxq; thisval];
                        
                        %Jfxq = [Jfxq; zeros(1,size(Jfxq,2))];
                        %Jfxq(end,(structure.state_order+1)*q) = 1;
                        
                    end
                    
                    % include Jacobian of controls, if present
                    
                    if structure.num_controls>0
                        
                        Jfuq = tmpJu{end}(1+(structure.state_order+1)*(q-1):(structure.state_order+1)*q,:);
                        %Jfuq(1,:) = zeros(1,size(Jfuq,2));
                        
                        if (structure.imposed_final_states(q)==1)
                            
                            Jfuq = [Jfuq; zeros(1,size(Jfuq,2))];
                            
                        end
                        
                        
                    end
                    
                    if structure.num_controls>0
                        
                        J_temp = [Jfxq Jfuq];
                        
                    else
                        
                        J_temp = Jfxq;
                        
                    end
                    
                    
                    endy = starty + size(J_temp,1)-1;
                    endy_aug = starty_aug +size(J_temp,1)-1;
                    endx_aug = startx_aug +size(J_temp,2)-1;
                    
                    J(starty_aug:endy_aug,startx_aug:endx_aug) = J_temp;
                    
                    % apply matching conditions to previous element
                    thisval = prevval(q,:);
                    J(starty_aug,(structure.num_elems-2)*(structure.num_eqs*(structure.state_order+1)+structure.num_controls*(structure.control_order+1))+1:(structure.num_elems-2)*(structure.num_eqs*(structure.state_order+1)+structure.num_controls*(structure.control_order+1))+length(thisval)) = thisval;
                    
                    starty = endy+1;
                    starty_aug = endy_aug+1;
                    
                end
                
            else
                
                startx = 1;
                starty = 1;
                
                % first element
                
                Jtmp = structure.M(1:(structure.test_order+1)*structure.num_eqs,1:(structure.state_order+1)*structure.num_eqs) + tmpJx{1};
                
                % add Jacobian of controls, if present
                if structure.num_controls>0
                    
                    Jtmp = [Jtmp tmpJu{1}];
                    
                end
                
                endx = startx+size(Jtmp,2)-1;
                endy = starty+size(Jtmp,1)-1;
                
                J(starty:endy,startx:endx) = Jtmp;
                
                starty = endy+1;
                startx = endx+1;
                
                % 2nd to last element
                
                ids = 1:size(tmpJx{1},1);
                local_rows_to_ignore = 1:structure.test_order+1:size(tmpJx{1},1);
                ids(local_rows_to_ignore) = [];
                
                startx_match = 0;
                
                starty_match = 0;
                
                for i = 2:structure.num_elems
                    
                    startx_match = startx_match + (structure.state_order+1)*structure.num_eqs+1*(i==2);
                    endx_match = startx_match+(structure.state_order+1)*structure.num_eqs-1;
                    
                    thisM = structure.M(1+(structure.test_order+1)*structure.num_eqs+structure.test_order*structure.num_eqs*(i-2):(structure.test_order+1)*structure.num_eqs+structure.test_order*structure.num_eqs*(i-1),1+(structure.state_order+1)*(i-1)*structure.num_eqs:(structure.state_order+1)*i*structure.num_eqs);
                    
                    Jtmp = tmpJx{i}(ids,:);
                    
                    Jtmp = thisM + Jtmp;
                    
                    % add Jacobian of controls, if present
                    if structure.num_controls>0
                        
                        Jtmp = [Jtmp tmpJu{i}(ids,:)];
                        
                    end
                    
                    endx = startx+size(Jtmp,2)-1;
                    endy = starty+size(Jtmp,1)-1;
                    
                    J(starty:endy,startx:endx) = Jtmp;
                    
                    % add matching condition on previous element
                    
                    for q = 1:structure.num_eqs
                        
                        starty_match = starty_match + structure.test_order+1*(i==2);
                        
                        thisM = structure.M(starty_match,startx_match:endx_match);
                        Jtmp = tmpJx{i}(local_rows_to_ignore(q),:);
                        Jtmp = thisM + Jtmp;
                        
                        if structure.num_controls>0
                            
                            Jtmp = [Jtmp tmpJu{i}(local_rows_to_ignore(q),:)];
                            
                        end
                        
                        J(starty_match,startx:endx) = Jtmp;
                        
                    end
                    
                    starty = endy+1;
                    startx = endx+1;
                    
                end
                
                % apply final conditions to last element
                
                if structure.num_free_states>0
                    
                    starty_match = starty_match+1;
                    startx = endx +1;
                    valr = structure.test_basis{end}(1)';
                    
                    for q = 1:structure.num_eqs
                        
                        if structure.free_final_states(q)==1
                            
                            J(starty_match:starty_match+structure.test_order-1,startx) = -valr(2+(structure.test_order+1)*(q-1):(structure.test_order+1)*q,q);
                            startx = startx+(structure.free_final_states(q)==1);
                            
                        end
                        
                        starty_match = starty_match + structure.test_order;
                        
                    end
                    
                end
                
            end
            
        end
        
    end
    
else
    
    J = [];
    
end

end