function [F,J] = eval_path_constraints(structure,x,x_0,x_f,t_0,t_f,static,els,calc_jac)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%
% This function evaluates inequality path constraints and their Jacobian


%% Extraction of nodal values of states and controls

ee1 = els(:,1)*(t_f-t_0);
ee2 = els(:,2)*(t_f-t_0);

loc_state = zeros((structure.state_order+1)*structure.num_eqs,structure.num_elems); % preallocation, will be needed to more easily map between the input vector (mixed states and controls) and the states
loc_u = zeros((structure.control_order+1)*structure.num_controls,structure.num_elems); % preallocation, will be needed to more easily map between the input vector (mixed states and controls) and the controls

startx = 1;

% mapping mixed states and controls of input x vector into more easily
% managed states and controls

for i=1:structure.num_elems
    
    endx = startx+(structure.state_order+1)*structure.num_eqs-1;
    loc_state(:,i) = x(startx:endx);
    
    startx = endx+1*(structure.num_controls>0);
    endx = startx+(structure.control_order+1)*structure.num_controls-1*(structure.num_controls>0);
    
    loc_u(:,i) = x(startx:endx);
    
    startx = endx+1;
    
end

%% Extract xb

[~,~,x_b] = extract_solution(x,structure,x_f);

%% Get scaling factors

scales = structure.scales;

%% Computation of the constraints on collocation points

num_con = length(structure.c(structure.state_eval_nodes{1,1}*loc_state(:,1),structure.control_eval_nodes{1,1}*loc_u(:,1),(ee1(1)+ee2(1))/2+(ee2(1)-ee1(1))/2*structure.col_nodes(1),static,scales,structure.constants));

% Path constraints will be evaluated at the (all the )nodes, so they
% will be imposed strongly in the collocation points. NOTE: if the
% order or type of polynomials is different between states and
% controls, some of the state nodes and control nodes will be
% mismatched. Inequalities will be evaluated in ALL these points,so the
% number of inequalities depends on the number of elements and the
% number of unique nodes.

lf = (size(structure.state_eval_nodes,1)*size(structure.state_eval_nodes,2)+(structure.DFET*2))*num_con;
F = zeros(lf,1);
start =1;

% Eval constraints on initial condition (only for DFET for now...)

if structure.DFET==1
    
    endd = start+num_con-1;
    
    F(start:endd) = structure.c(x_0,structure.control_valsl*loc_u(:,1),t_0,static,scales,structure.constants);
            
    start = endd+1;

end

% Eval constraints on central nodes

for i = 1:structure.num_elems
    
    e1 = ee1(i);
    e2 = ee2(i);
    
    for j = 1:length(structure.col_nodes)
        
        endd = start+num_con-1;
        
        F(start:endd) = structure.c(structure.state_eval_nodes{i,j}*loc_state(:,i),structure.control_eval_nodes{i,j}*loc_u(:,i),(e2+e1)/2+(e2-e1)/2*structure.col_nodes(j),static,scales,structure.constants);
        
        start = endd+1;
        
    end
    
end

% add constraints also to final nodes, in case of DFET

if structure.DFET==1
    
    endd = start+num_con-1;
    
    F(start:endd) = structure.c(x_b,structure.control_valsr*loc_u(:,i),t_f,static,scales,structure.constants);
    
end

%% Computation of the Jacobian in the projected state space

J = zeros(length(F),length(x));

if calc_jac
    
%    if isempty(structure.dcx)
%        
%         %% FINITE DIFFERENCES
%         
%         % first order "forward" finite difference approach with fixed step
%         % of 1e-6. Would like to use backwards differences when variable is
%         % at upper bound, but would require to know the bounds here. Second
%         % order accurate finite differences could also be appealing,
%         % requiring to average forward and backward differences and
%         % eventually using ad adapted stencil for the boundaries, but
%         % probably the additional computational cost does not pay off.
%         
%         for i = 1:length(x)
%             
%             x_v = x;
%             x_v(i) = x_v(i)+0.000001;
%             F_temp = eval_path_constraints(structure,x_v,x_0,x_f,t_0,t_f,els,0);
%             J(:,i) = (F_temp-F)/0.000001;
%             
%         end
        
%    else
        
        %% ANALYTIC JACOBIAN
        
        ypos =1;
        startx = 1;
        dx = size(structure.state_eval_nodes{1,1},2)+size(structure.control_eval_nodes{1,1},2);
        endx = startx+dx-1;
       
        % derivativs of constraints evaluated at initial condition. Only
        % controls derivatives appear here, because derivs wrt free states
        % are performed outside (as they are static variables). Derivatives
        % wrt states are computed and then zeroed out to keep the size of
        % the contribution correct
        
        if structure.DFET==1
                        
            dcx = structure.dcx(x_0,structure.control_valsl*loc_u(:,1),t_0,static,scales,structure.constants)*structure.state_valsl;
            dcu = structure.dcu(x_0,structure.control_valsl*loc_u(:,1),t_0,static,scales,structure.constants)*structure.control_valsl;
            
            deriv = [dcx*0 dcu];
            
            J(ypos:ypos+num_con-1,startx:endx) = deriv;
            
            ypos = ypos+num_con;
            
        end
        
        % derivatives of constraints evaluated at intermediate points

        for i = 1:structure.num_elems
            
            endx = startx+dx-1;
            
            e1 = ee1(i);
            e2 = ee2(i);
           
            for j = 1:length(structure.col_nodes)
                                
                dcx = structure.dcx(structure.state_eval_nodes{i,j}*loc_state(:,i),structure.control_eval_nodes{i,j}*loc_u(:,i),(e2+e1)/2+(e2-e1)/2*structure.col_nodes(j),static,scales,structure.constants)*structure.state_eval_nodes{i,j};
                dcu = structure.dcu(structure.state_eval_nodes{i,j}*loc_state(:,i),structure.control_eval_nodes{i,j}*loc_u(:,i),(e2+e1)/2+(e2-e1)/2*structure.col_nodes(j),static,scales,structure.constants)*structure.control_eval_nodes{i,j};
                
                J(ypos:ypos+num_con-1,startx:endx) = [dcx dcu];
                
                ypos = ypos+num_con;
                
            end
            
            startx = endx+1;
            
        end

        % derivativs wrt final nodes, in case of DFET (remember, in that
        % inequality also controls must be added)
        
        if structure.DFET==1
                        
            dcx = structure.dcx(x_b,structure.control_valsr*loc_u(:,i),t_f,static,scales,structure.constants);%*structure.state_eval_nodes{i,j};
            dcu = structure.dcu(x_b,structure.control_valsr*loc_u(:,i),t_f,static,scales,structure.constants)*structure.control_valsr;

            % only keep derivatives wrt free final values
            
            dcx = dcx(:,structure.free_final_states);
            
            deriv = [dcu dcx];
            
            len = size(deriv,2);

            J(ypos:ypos+num_con-1,end-len+1:end) = deriv;
            
        end
        
%    end
    
else
    
    J = [];
    
end

J = sparse(J);

end