function [val,x_sol] = max_ener_MACS_MOO(x_in,lb,ub,f,structure,x_0,x_f,dfx,dfu,dft,fminconoptions)

[x_sol,~,exitflag,~,~] = fmincon(@(x) 1,x_in(2:end),[],[],[],[],lb(2:end),ub(2:end),@(x) dynamics2(f,structure,x,x_0,x_f,x_in(1),dfx,dfu),fminconoptions);

x_sol = [x_in(1) x_sol];

if exitflag==1
    
    t_f = x_sol(1);
    t_0 = 0;
    
    g = @(x,u,t) [t 0; -(x(2)^2+x(4)^2)/2+1/x(1) 0];    
    %dgdx = @(x,u,t) [0 0];
    %dgdu = @(x,u,t) [0 0];
    %dgdt = @(x,u,t) [1 0];
    weights = [1 0; 1 0];
    [x_out,u,x_b] = extract_solution(x_sol(2:end),structure,x_f);
    
    % finite differences
    %[~,F0,~,J0] = dynamics2(f,structure,x_sol,x_0,x_f,dfx,dfu);
    
    %J0 = J0(2:end,:)';
  
    %[dFdx,dFdu,dFdxb] = extract_derivatives(J0,structure);
    
    %epsi = 1e-9;
    %xtrial = x_sol;
    %xtrial(1) = xtrial(1)+epsi;
    %[~,diff] = dynamics2(f,structure,xtrial,x_0,x_f,dfx,dfu);
        
    %dFdtf = (diff-F0)/epsi;

    %gradu = [1;pinv(dFdu)*dFdtf]; %needs to be "rearranged" to be ordered as xsol
    
    %gradx = [1;pinv(dFdx)*dFdtf];
    
    %grad = -[1;pinv(J0)*dFdtf]';
    
    %grad = [-1 zeros(1,length(x_sol)-1)];
    
    %dg = 1;
    
    %dudtf = pinv(dFdu)*dFdtf;
    %dxfdtf = pinv(dFdxb)*dFdtf;
    
    %sens = eval_sensitivity(f,structure,x_sol,x_0,x_f,t_0,t_f,dfx,dfu,dft);        
     
    val = eval_cost_functions(g,weights,x_out,u,x_b,[t_0 t_f],structure,0,[],[])';
    
    %dval = [1 dxfdtf(2)];
    
else
    
    %grad = nan(1,length(x_sol));
    
    val = [1e9 1e9];
    %dval = [0];% 0]';
end

end

function [c,ceq,Jc,Jceq] = dynamics2(f,structure,x_in,x_0,x_f,t_f,dfx,dfu)

c = [];
Jc = [];

t_0 = 0;
%t_f = x_in(1);
x_guess = x_in;%(2:end);

[ceq,Jceq] = eval_constraints(f,structure,x_guess,x_0,x_f,t_0,t_f,1,dfx,dfu);

Jceq = Jceq';
%f2 = eval_constraints(f,structure,x_guess,x_0,x_f,t_0,t_f+1e-6,0,[],[]);
%dfdt = (f2-ceq)'/1e-6;
%Jceq = [dfdt;Jceq];

end