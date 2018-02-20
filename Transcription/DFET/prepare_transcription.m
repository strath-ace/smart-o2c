function out = prepare_transcription(num_eqs,num_controls,num_elems,state_order,control_order,integr_order,DFET,state_distrib,control_distrib)
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%-----------Copyright (C) 2016 University of Strathclyde-------------
%
%
%
% Transcribes problem
% BLA BLA BLA later on...
% Mi aspetto f(x,u,t) direttamente, quindi num_eqs e num_controls sono dati
% direttamente dall'utente per evitare casino

%% Inputs check (TOFINISH)

if (DFET ~=0) && (DFET~=1)
    
    error('DFET must be 0 (continuous) or 1 (discontinuous)')
    
end

if (~strcmp(state_distrib,'Legendre'))&&(~strcmp(state_distrib,'Cheby')&&(~strcmp(state_distrib,'Lobatto')))

    error('state_distrib must be Cheby, Legendre or Lobatto')
    
end

if (~strcmp(control_distrib,'Legendre'))&&(~strcmp(control_distrib,'Cheby')&&(~strcmp(control_distrib,'Lobatto')))

    error('control_distrib must be Cheby, Legendre or Lobatto')
    
end

%% helper variables

test_order = state_order+(DFET==1);

%% Construction of vector of nodes (assume equispaced at this stage, but arbitrary spacing is allowed)
% options will contain a flag to change element type (i.e. linear or Chebychev distribution of nodes inside each element)

t_n = linspace(0,1,num_elems+1);    % this are the "external" nodes of each element

[els, in_nodes_state,s_nodes] = make_elements(t_n,state_order,state_distrib);
[~, in_nodes_control,c_nodes] = make_elements(t_n,control_order,control_distrib);
[~, ~,t_nodes] = make_elements(t_n,test_order,'Lobatto');   %always with lobatto because it simplifies a lot DFET, not needed for continuous although it helps even there

%% Computation of weights and nodes for Gauss quadrature
% options will contain a flag to change integration scheme (Gauss-Legendre or Gauss Lobatto)

[integr_nodes,int_weights] = lgwt(ceil((integr_order+1)/2),-1,1);   % Gauss-Legendre nodes
%[int_nodes,int_weights] = lgwt(2*(order+1),-1,1);  % Gauss-Lobatto nodes

%% Construction of basis directly through Lagrange Polynomials, derivatives

[PC,DPC,state_basis,state_eval] = make_basis2(state_order,num_elems,num_eqs,s_nodes,integr_nodes);
[PCU,~,control_basis,control_eval] = make_basis2(control_order,num_elems,num_controls,c_nodes,integr_nodes);
[state_eval_nodes,control_eval_nodes,col_nodes] = eval_poly_on_nodes(num_elems,state_basis,s_nodes,control_basis,c_nodes);


if (DFET==1)
    % Basis for test function!!!! (in DGFET is one order higher than nonDGFET)
    [PCT,DPCT,test_basis,test_eval] = make_basis2(test_order,num_elems,num_eqs,t_nodes,integr_nodes);
    
else
   
    PCT = PC;
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
out.s_nodes = s_nodes;
out.c_nodes = c_nodes;
out.test_basis = test_basis;
out.test_eval = test_eval;
out.test_order = test_order;
out.state_eval_nodes = state_eval_nodes;
out.control_eval_nodes = control_eval_nodes;
out.col_nodes = col_nodes;

end
