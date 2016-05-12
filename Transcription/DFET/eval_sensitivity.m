function sens = eval_sensitivity(f,structure,x,x_0,x_f,t_0,t_f,dfx,dfu,dft)

% sens is a structure containing the sensitivities of a CONVERGED
% (FEASIBLE) solution of the system of equations (i.e. both the "x" and the
% "u") wrt all parameters (x_0, x_f, t_0, t_f and "qu", i.e. coefficients 
% of all the controls)

loc_state = zeros((structure.state_order+1)*structure.num_eqs,structure.num_elems);
loc_u = zeros((structure.control_order+1)*structure.num_controls,structure.num_elems);

startx = 1;

% extracts loc_state and loc_u from x, making following computations much
% clearer. Shouldn't influence the computation of derivatives. It's done
% once for all the sub matrices

for i=1:structure.num_elems
        
    endx = startx+(structure.state_order+1)*structure.num_eqs-1;
    loc_state(:,i) = x(startx:endx);
    
    startx = endx+1*(structure.num_controls>0);
    endx = startx+(structure.control_order+1)*structure.num_controls-1*(structure.num_controls>0);
    
    loc_u(:,i) = x(startx:endx);
    
    startx = endx+1;
    
end

% these contain the derivatives of ALL states, I need to extract only the
% final states... should be immediate, but it doesn't seem so...

sens.dFdt0 = eval_sensitivity_wrt_t0(f,structure,t_0,t_f,loc_state,loc_u,dft);
sens.dFdtf = eval_sensitivity_wrt_tf(f,structure,t_0,t_f,loc_state,loc_u,dft);


end