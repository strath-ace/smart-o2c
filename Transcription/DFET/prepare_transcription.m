function out = prepare_transcription(num_eqs,num_controls,num_elems,state_order,control_order,DFET,state_distrib,control_distrib,test_distrib,integr_type)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%
% Creates elements, nodes, mass matrix, basis functions, test function, and
% all the underlying machinery of DFET

%% Inputs check (TOFINISH)

if (DFET ~=0) && (DFET~=1)
    
    error('DFET must be 0 (continuous) or 1 (discontinuous)')
    
end

if (~strcmp(state_distrib,'Legendre'))&&(~strcmp(state_distrib,'Bernstein')&&(~strcmp(state_distrib,'Lobatto')))

    error('state_distrib must be Bernstein, Legendre or Lobatto')
    
end

if (~strcmp(control_distrib,'Legendre'))&&(~strcmp(control_distrib,'Bernstein')&&(~strcmp(control_distrib,'Lobatto')))

    error('control_distrib must be Bernstein, Legendre or Lobatto')
    
end

if (~strcmp(test_distrib,'Legendre'))&&(~strcmp(test_distrib,'Bernstein')&&(~strcmp(test_distrib,'Lobatto')))

    error('test_distrib must be Bernstein, Legendre or Lobatto')
    
end

if (~strcmp(integr_type,'Legendre'))&&(~strcmp(integr_type,'Lobatto'))

    error('integr_type must be Legendre or Lobatto')
    
end

if num_elems==2
   
    error('There is a bug when using exactly 2 elements in the computation of the Jacobian of the dynamics. Will have to fix it. Meanwhile use more elements');
    
end

%% helper variables

test_order = state_order+(DFET==1);

%% Construction of vector of nodes (assume equispaced at this stage, but arbitrary spacing is allowed)
% options will contain a flag to change element type (i.e. linear or Chebychev distribution of nodes inside each element)

t_n = linspace(0,1,num_elems+1);    % this are the "external" nodes of each element

[els, in_nodes_state,s_nodes] = make_elements(t_n,state_order,state_distrib);
[~, in_nodes_control,c_nodes] = make_elements(t_n,control_order,control_distrib);
[~, ~,t_nodes] = make_elements(t_n,test_order,test_distrib);   

%% Computation of weights and nodes for Gauss quadrature

switch (integr_type)
    
    case 'Legendre'
        
        [integr_nodes,int_weights] = lgwt(max(test_order,control_order),-1,1); % Gauss-Lobatto nodes
        
    case 'Lobatto'
        
        [integr_nodes,int_weights] = lglnodes(max(test_order-1,control_order));  % Gauss-Lobatto nodes
        
%     case 'Radau'
%     
%         [integr_nodes,int_weights] = lgrnodes(max(test_order,control_order));  % Gauss-Radau nodes

end

integr_nodes = fliplr(integr_nodes(:)');

%% Construction of basis directly through Lagrange Polynomials, derivatives

[PC,DPC,state_basis,state_eval] = make_basis2(state_order,num_elems,num_eqs,s_nodes,integr_nodes,state_distrib);
[PCU,~,control_basis,control_eval] = make_basis2(control_order,num_elems,num_controls,c_nodes,integr_nodes,control_distrib);
[state_eval_nodes,control_eval_nodes,col_nodes] = eval_poly_on_nodes(num_elems,state_basis,control_basis,s_nodes);

if (DFET==1)
    % Basis for test function!!!! (in DGFET is one order higher than nonDGFET)
    [~,DPCT,test_basis,test_eval] = make_basis2(test_order,num_elems,num_eqs,t_nodes,integr_nodes,test_distrib);
    
else
   
    %PCT = PC;
    DPCT = DPC;
    test_basis = state_basis;
    test_eval = state_eval;
    
end

%% Construction of Mass Matrix M

%pos = ones(state_order+test_order+1,1);
pos = ones(max(size(DPCT,1),1)+size(PC,1),1);
expon = (length(pos)-1:-1:0)';
neg = (-1).^expon;

M = make_mass_matrix(num_elems,state_order,test_order,num_eqs,PC,DPCT,pos,neg,DFET);

%% Assign outputs

out.M = M;
out.state_eval = state_eval;
out.control_eval = control_eval;
out.uniform_els = els;
out.uniform_in_nodes_state = in_nodes_state;
out.uniform_in_nodes_control = in_nodes_control;
out.integr_nodes = integr_nodes;
out.integr_weights = int_weights;
out.PC = PC;
out.DPC = DPC;
out.state_basis = state_basis;
out.state_eval = state_eval;
out.PCU = PCU;
out.control_basis = control_basis;
out.control_eval = control_eval;
out.DFET = DFET;
out.state_order = state_order;
out.control_order = control_order;
out.num_elems = num_elems;
out.num_eqs = num_eqs;
out.num_controls = num_controls;
out.state_distrib = state_distrib;
out.control_distrib = control_distrib;
out.test_distrib = test_distrib;
out.s_nodes = s_nodes;
out.c_nodes = c_nodes;
out.test_basis = test_basis;
out.test_eval = test_eval;
out.test_order = test_order;
out.state_eval_nodes = state_eval_nodes;
out.control_eval_nodes = control_eval_nodes;
out.col_nodes = col_nodes;

%% precalculating often used values and saving them as sparse matrices
% left value of test functions on first element
tmp = test_basis{1}(-1)';
tmp(abs(tmp)<1e-9) = 0;         % chopping numerical roundoff errors
tmp = sparse(tmp);
out.basis_valsl1 = tmp;

% left value of test functions on last element
tmp = test_basis{end}(-1)';
tmp(abs(tmp)<1e-9) = 0;         % chopping numerical roundoff errors
tmp = sparse(tmp);
out.basis_valsl =  tmp;

% right value of test functions on last element
tmp = test_basis{end}(1)';
tmp(abs(tmp)<1e-9) = 0;         % chopping numerical roundoff errors
tmp = sparse(tmp);
out.basis_valsr = tmp;

% left value of state basis functions on first element
tmp = state_basis{1}(-1);
tmp(abs(tmp)<1e-9) = 0;         % chopping numerical roundoff errors
tmp = sparse(tmp);
out.state_valsl = tmp;

% right value of state basis functions on last element
tmp = state_basis{end}(1);
tmp(abs(tmp)<1e-9) = 0;         % chopping numerical roundoff errors
tmp = sparse(tmp);
out.state_valsr = tmp;

% right value of controls basis functions on last element
tmp = control_basis{end}(1);
tmp(abs(tmp)<1e-9) = 0;         % chopping numerical roundoff errors
tmp = sparse(tmp);
out.control_valsr = tmp;

% left value of controls basis functions on first element
tmp = control_basis{1}(-1);
tmp(abs(tmp)<1e-9) = 0;         % chopping numerical roundoff errors
tmp = sparse(tmp);
out.control_valsl = tmp;

end