function [vals,grad] = eval_cost_functions2(g,weights,x_in,x_0,x_f,time_interval,static,els,structure,compute_grad,dgu0,dgxf,dguf,dgxi,dgui)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%
% Evaluates generic cost functions in the Bolza form

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
%
% grad gives the gradients of all objective functions wrt NODAL
% states and controls. They are ordered the same way as the objective
% functions themselves, thus first column of grad contains gradients of
% g(t_f) while second column contains gradients of int(g)

%% If no objective function is defined, just return 0

if isempty(g)
   
    vals = 0;
    grad = zeros(length(x_in),1);
    
   return 
   
end

%% Get scaling factors

scales = structure.scales;

%% Check input

[num_cost_funcs,n_cols] = size(g(zeros(structure.num_eqs,1),zeros(structure.num_controls,1),0,zeros(structure.num_eqs,1),zeros(structure.num_controls,1),0,zeros(structure.num_eqs,1),zeros(structure.num_controls,1),0,static,scales,structure.constants));    % number of functions

if n_cols~=2
    
    error('Wrong format for g. Must be a function handle of (x,u,t), with exactly 2 columns, and as many rows as the number of cost functions')
    
end

if size(weights,1)~=num_cost_funcs
    
    error('Weights vector must have same number of columns as g')
    
end

if size(weights,2)~=2
    
    error('Weights must have exactly 2 columns')
    
end

if any(all(weights==0,2))
    
    error('At least one row of the weights vector is identically 0')
    
end

if length(time_interval)~=2
    
    error('interval must be a vector containing t_0 and t_f only');
    
end

if (compute_grad~=0) && (compute_grad~=1)
    
    error('compute_grad can be either 0 or 1');
    
end

% if compute_grad==1
%
%     if ~all([isempty(dgxf),isempty(dguf),isempty(dgxi),isempty(dgui)])
%
%         error('Either give both gradients of c(t_f) and int(c) wrt states and controls or give none');
%
%     end
%
%     if ~isempty(dgx)
%
%         error('Not implemented, yet..')
%
% %         if size(dgx(zeros(structure.num_eqs,1),zeros(structure.num_controls,1),0),1)~=num_cost_funcs
% %
% %             error('dgx must have same number of rows as g')
% %
% %         end
% %
% %         if size(dgu(zeros(structure.num_eqs,1),zeros(structure.num_controls,1),0),1)~=num_cost_funcs
% %
% %             error('dgu must have same number of rows as g')
% %
% %         end
%
%     end
%
% end

t_0 = time_interval(1);%min(time_interval);
t_f = time_interval(2);%max(time_interval);
% 
% if t_0==t_f
%     
%     error('t_0 and t_f must be different!');
%     
% end

%% Extract x,u,xb

[x,u,x_b] = extract_solution(x_in,structure,x_f);

u_0 = structure.control_valsl*u(:,1);
u_f = structure.control_valsr*u(:,end);

%% Compute values

vals = [];

if any(weights(:,1)>0)
    
    % evaluate "boundary value functions", and multiply them by their weights
    
    if structure.DFET==0
        
        qq = g(structure.state_valsr*x(:,end),structure.control_basis{end}(1)*u(:,end),t_f,static,scales,structure.constants);
        
    else
        
        % the variables for the integral part need to be passed, even if
        % it's not evaluated. We pass arbitrarily the final values
        
        qq = g(x_0,u_0,t_0,x_b,u_f,t_f,x_b,u_f,t_f,static,scales,structure.constants);
        
    end
    
    vals = qq(:,1).*weights(:,1);
    
end

if any(weights(:,2)>0)
    
    if isempty(vals)
        
        vals = 0;
        
    end
    
    % compute integrals elementwise (can be made PARALLEL, although benefits depend widely upon the actual amount of work done by each node)
    
    ee1 = els(:,1)*(t_f-t_0)*structure.scales.tscale;
    ee2 = els(:,2)*(t_f-t_0)*structure.scales.tscale;
    
    for i = 1:structure.num_elems
        
        e1 = ee1(i);
        e2 = ee2(i);
        
        for j = 1:length(structure.integr_nodes)
            
            qq = structure.integr_weights(j)*g(x_0,u_0,t_0,x_b,u_f,t_f,structure.state_eval{i,j}*x(:,i),structure.control_eval{i,j}*u(:,i),(e2+e1)/2+(e2-e1)/2*structure.integr_nodes(j),static,scales,structure.constants).*(e2-e1)/2;
            
            vals = vals + qq(:,2).*weights(:,2);
            
        end
        
    end
    
end

%% Compute gradients

grad = zeros(length(x_in),num_cost_funcs);

if compute_grad==1
    
%    if isempty(dgxf)
%        
%         %% FINITE DIFFERENCES
%         
%         grad = zeros(num_cost_funcs,length(x_in));
%         
%         % second order "centered" finite difference approach with fixed
%         % step h
%         h=1e-9;
%         
%         for i = 1:length(x_in)
%             
%             x_v = x_in;
%             x_v(i) = x_v(i)+h;
%             g_temp = eval_cost_functions2(g,weights,x_v,x_f,time_interval,els,structure,0,[],[],[],[]);
%             
%             x_v = x_in;
%             x_v(i) = x_v(i)-h;
%             g_temp2 = eval_cost_functions2(g,weights,x_v,x_f,time_interval,els,structure,0,[],[],[],[]);
%             
%             
%             grad(:,i) = (g_temp-g_temp2)/(2*h);
%             
%         end
        
%    else
        
        if any(weights(:,1)>0)
            
            % evaluate "final state" gradients, and multiply by their weights
            
            if structure.DFET==0
                
                error('Gradients of objectives not implemented for bi-continuous DFET, yet...')
                
                %             % Derivatives wrt u(t_f) instead are computed by projection
                %             % onto the control basis of the final element, and stored
                %             % accordingly.
                %             % TODO: these gradients practically never appear, so computing
                %             % them when they are always 0 is a waste. Need to find an
                %             % effective approach to this. Flag?
                %
                %             graduf = (dguf(x_b,structure.control_valsr*u(:,end),t_f)*structure.control_valsr.*weights(:,1))';
                %
                %             % place graduf in the right position
                %
                %             xend = length(x_in);
                %             xstart = xend-size(graduf)+1;
                %
                %             grad(xstart:xend,:) = graduf;
                
            else
                
                % Derivatives wrt u(t_0)

                gradu0 = (dgu0(x_0,u_0,t_0,x_b,u_f,t_f,x_b,u_f,t_f,static,scales,structure.constants)*structure.control_valsl.*weights(:,1))';
                
                % place graduf in the right position
                                
                xstart = size(structure.state_eval_nodes{1,1},2)+1;
                xend = xstart+size(gradu0,1)-1;
                
                grad(xstart:xend,:) = gradu0;
                
                % Derivatives wrt x_f
                
                if structure.num_free_states>0
                    
                    % extract id of free final states
                    
                    vect = 1:length(x_f);
                    ids = vect(structure.free_final_states);
                    
                    % eval gradient function (returns gradient wrt all final
                    % states) and multiply it to the weights vector
                    % note: gradxf has the same ordering as the Jacobian of the
                    % dynamics (objective functions are sorted by row,
                    % variables by column), so we need to transpose it to match
                    % the ordering required by fmincon
                    
                    gradxf = (dgxf(x_0,u_0,t_0,x_b,u_f,t_f,x_b,u_f,t_f,static,scales,structure.constants).*weights(:,1))';
                    
                    % extracting only compotents wrt free final states
                    
                    gradxf = gradxf(ids,:);
                    
                    % place those components at their place, i.e. at the very
                    % end of the solution vector
                    
                    grad(end-structure.num_free_states+1:end,:) = gradxf;
                    
                end
                
                % Derivatives wrt u(t_f) 

                graduf = (dguf(x_0,u_0,t_0,x_b,u_f,t_f,x_b,u_f,t_f,static,scales,structure.constants)*structure.control_valsr.*weights(:,1))';
                
                % place graduf in the right position
                
                xend = length(x_in)-structure.num_free_states;
                xstart = xend-size(graduf)+1;
                
                grad(xstart:xend,:) = graduf;
                
            end
            
            
        end
        
        if any(weights(:,2)>0)
            
            % gradients (wrt states and controls) of the integral of objective
            % functions is the integral of the gradients (wrt states and
            % controls) of the objective functions
            
            if structure.DFET==0
                
                error('Gradients of objectives not implemented for bi-continuous DFET, yet...')
                
            else
                
                % for each element
                
                startx = 1;
                
                for i = 1:structure.num_elems
                    
                    e1 = ee1(i);
                    e2 = ee2(i);
                    
                    gx = sparse(length(structure.state_eval{1,1}),num_cost_funcs);
                    gu = sparse(length(structure.control_eval{1,1}),num_cost_funcs);
                    
                    % perform integration inside this element (as before, those
                    % values have to be transposed to be consistent)
                    
                    for j = 1:length(structure.integr_nodes)
                        
                        gx = gx + (structure.integr_weights(j)*dgxi(x_0,u_0,t_0,x_b,u_f,t_f,structure.state_eval{i,j}*x(:,i),structure.control_eval{i,j}*u(:,i),(e2+e1)/2+(e2-e1)/2*structure.integr_nodes(j),static,scales,structure.constants)*structure.state_eval{i,j}.*(e2-e1)/2)';
                        gu = gu + (structure.integr_weights(j)*dgui(x_0,u_0,t_0,x_b,u_f,t_f,structure.state_eval{i,j}*x(:,i),structure.control_eval{i,j}*u(:,i),(e2+e1)/2+(e2-e1)/2*structure.integr_nodes(j),static,scales,structure.constants)*structure.control_eval{i,j}.*(e2-e1)/2)';
                        
                    end
                    
                    % place integrals at the right position in grad
                    
                    endx = startx+size(gx,1)-1;
                    grad(startx:endx,:) = grad(startx:endx,:)+gx.*structure.weights(:,2)';
                    
                    startx = endx+1;
                    endx = startx+size(gu,1)-1;
                    grad(startx:endx,:) = grad(startx:endx,:)+gu.*structure.weights(:,2)';
                    
                    startx = endx+1;
                    
                end
                
            end
            
        end
        
    %end
    
end