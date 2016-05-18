function [vals,gradqu] = eval_cost_functions2(g,weights,x,u,x_b,time_interval,structure,compute_grad,dgxf,dguf,dgxi,dgui)

% g is a 2 columns vector containing function handles. left entries are
% for functions for which we are interested in evaluating them at final
% time, right entry is for functions which we are interested in evaluating
% their INTEGRAL between t0 and tf
% weights is a vector (same size as g), to compute the weighted sum between
% the two values
%
% EXAMPLE: Bolza problem 1
%
% g = [g_1 g_2], weights = [1/2 1/2].
%
% This means that the final val= g_1(tf)/2 + int(g_2,t0,tf)/2
%
% x and u vectors must be already reordered element by element


%% Check input

[num_cost_funcs,n_cols] = size(g(zeros(structure.num_eqs,1),zeros(structure.num_controls,1),0));

if n_cols~=2
    
    error('Wrong format for g. Must be a function handle of (x,u,t), with exactly 2 columns, and as many rows as the number of cost functions')
    
end

if size(weights,1)~=num_cost_funcs
    
    error('Weights vector must have same number of columns as g')
    
end

if size(weights,2)~=2
    
    error('Weights must have exactly 2 columns')
    
end

if length(time_interval)~=2
    
    error('interval must be a vector containing t_0 and t_f only');
    
end

if (compute_grad~=0) && (compute_grad~=1)
    
    error('compute_grad can be either 0 or 1');
    
end

if compute_grad==1
    
    if isempty(dgxf)|| isempty(dguf) || isempty(dgxi)|| isempty(dgui)
        
        error('All derivatives of g(t_f) and int(g) wrt states and controls must be supplied');
        
    end
    
    if size(dgxf(zeros(structure.num_eqs,1),zeros(structure.num_controls,1),0),1)~=num_cost_funcs
        
        error('dgxf must have same number of rows as g')
        
    end
    
    if size(dgxf(zeros(structure.num_eqs,1),zeros(structure.num_controls,1)),2)~=structure.num_eqs
        
        error('dgxf must as many columns as the number of state variables (derivatives wrt all states)');
        
    end
    
    if size(dguf(zeros(structure.num_eqs,1),zeros(structure.num_controls,1),0),1)~=num_cost_funcs
        
        error('dguf must have same number of rows as g')
        
    end
    
    if size(dguf(zeros(structure.num_eqs,1),zeros(structure.num_controls,1),0),2)~=structure.num_controls
        
        error('dguf must as many columns as the number of control variables (derivatives wrt all controls)');
        
    end
    
    if size(dgxi(zeros(structure.num_eqs,1),zeros(structure.num_controls,1),0),1)~=num_cost_funcs
        
        error('dgxi must have same number of rows as g')
        
    end
    
    if size(dgxi(zeros(structure.num_eqs,1),zeros(structure.num_controls,1),0),2)~=structure.num_eqs
        
        error('dgxi must as many columns as the number of state variables (derivatives wrt all states)');
        
    end
    
    if size(dgui(zeros(structure.num_eqs,1),zeros(structure.num_controls,1),0),1)~=num_cost_funcs
        
        error('dgui must have same number of rows as g')
        
    end
    
    if size(dgui(zeros(structure.num_eqs,1),zeros(structure.num_controls,1),0),2)~=structure.num_controls
        
        error('dgui must as many columns as the number of control variables (derivatives wrt all controls)');
        
    end
    
end

t_0 = min(time_interval);
t_f = max(time_interval);

if t_0==t_f
    
    error('t_0 and t_f must be different!');
    
end

%% Compute values

vals = zeros(num_cost_funcs,1);

% evaluate "final state functions", and multiply them by their weights
if structure.DFET==0
    
    qq = g(structure.state_basis{end}(1)*x(:,end),structure.control_basis{end}(1)*u(:,end),t_f);
    
else
    
    qq = g(x_b,structure.control_basis{end}(1)*u(:,end),t_f);
    
end

vals = qq(:,1).*weights(:,1);

% compute integrals elementwise (can be made PARALLEL, although benefits depend widely upon the actual amount of work done by each node)

ee1 = structure.els(:,1);
ee2 = structure.els(:,2);

for i = 1:structure.num_elems
    
    e1 = ee1(i);
    e2 = ee2(i);
    
    % maybe it's better to avoid reshape, saves a for loop afterwards for
    % the reconstruction of complete f
    
    for j = 1:length(structure.integr_nodes)
        
        qq = structure.integr_weights(j)*g(structure.state_eval{i,j}*x(:,i),structure.control_eval{i,j}*u(:,i),(e2+e1)/2+(e2-e1)/2*structure.integr_nodes(j)).*(e2-e1)/2;
        
        vals = vals + qq(:,2).*weights(:,2);
        
    end
    
end

%% Compute gradients

if compute_grad==1
    
    gradx = zeros(1,(structure.state_order+1)*structure.num_eqs*structure.num_elems);
    gradu = zeros(1,(structure.control_order+1)*structure.num_controls*structure.num_elems);

    % evaluate "final state" gradients, and multiply by their weights
    
    if structure.DFET==0
                
        for i = 1:structure.num_elems
            
            qqxf = dgxf(structure.state_basis{end}(1)*x(:,end),structure.control_basis{end}(1)*u(:,end),t_f)*structure.state_basis{end}(1);
            qquf = dguf(structure.state_basis{end}(1)*x(:,end),structure.control_basis{end}(1)*u(:,end),t_f)*structure.control_basis{end}(1);
            
            gradx(1+(structure.state_order+1)*structure.num_eqs*(i-1):(structure.state_order+1)*structure.num_eqs*i) = repmat(weights(:,1),1,size(qqxf,2)).*qqxf;
            gradu(1+(structure.control_order+1)*structure.num_controls*(i-1):(structure.control_order+1)*structure.num_controls*i) = repmat(weights(:,1),1,size(qquf,2)).*qquf;
            
        end
        
    else
        
        for i = 1:structure.num_elems                        
            
            qqxf = dgxf(x_b,structure.control_basis{end}(1)*u(:,end),t_f)*structure.state_basis{end}(1);
            qquf = dguf(x_b,structure.control_basis{end}(1)*u(:,end),t_f)*structure.control_basis{end}(1);
            
            gradx(1+(structure.state_order+1)*structure.num_eqs*(i-1):(structure.state_order+1)*structure.num_eqs*i) = repmat(weights(:,1),1,size(qqxf,2)).*qqxf;
            gradu(1+(structure.control_order+1)*structure.num_controls*(i-1):(structure.control_order+1)*structure.num_controls*i) = repmat(weights(:,1),1,size(qquf,2)).*qquf;
            
        end
        
    end
    
    % evaluate "integral" gradients, and multiply by their weights
    
    for i = 1:structure.num_elems
        
        e1 = ee1(i);
        e2 = ee2(i);
        
        % maybe it's better to avoid reshape, saves a for loop afterwards for
        % the reconstruction of complete f
        
        for j = 1:length(structure.integr_nodes)
            
            qqxi = structure.integr_weights(j)*dgxi(structure.state_eval{i,j}*x(:,i),structure.control_eval{i,j}*u(:,i),(e2+e1)/2+(e2-e1)/2*structure.integr_nodes(j))*structure.state_eval{i,j}.*(e2-e1)/2;
            qqui = structure.integr_weights(j)*dgui(structure.state_eval{i,j}*x(:,i),structure.control_eval{i,j}*u(:,i),(e2+e1)/2+(e2-e1)/2*structure.integr_nodes(j))*structure.control_eval{i,j}.*(e2-e1)/2;
            
            gradx(1+(structure.state_order+1)*structure.num_eqs*(i-1):(structure.state_order+1)*structure.num_eqs*i) = gradx(1+(structure.state_order+1)*structure.num_eqs*(i-1):(structure.state_order+1)*structure.num_eqs*i) + qqxi.*repmat(weights(:,2),1,size(qqxi,2));
            gradu(1+(structure.control_order+1)*structure.num_controls*(i-1):(structure.control_order+1)*structure.num_controls*i) = gradu(1+(structure.control_order+1)*structure.num_controls*(i-1):(structure.control_order+1)*structure.num_controls*i) + qqui.*repmat(weights(:,2),1,size(qqui,2));
            
        end
        
    end
    
end

