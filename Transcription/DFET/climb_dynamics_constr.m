function [c,ceq,Jc,Jceq] = climb_dynamics_constr (f,g,structure,x_in,x_0,x_f,t_0,dfx,dfu)

% min time formulation
t_f = x_in(1);
x = x_in(2:end);

% max vel formulation
%t_f = 200;
%x = x_in;

if ~isempty(g)

    [c,Jc] = eval_path_constraints(g,structure,x,x_0,x_f,t_0,t_f,structure.uniform_els,1,1,[],[]);
   
    Jc = Jc';

    % min time formulation

    f2 = eval_path_constraints(g,structure,x,x_0,x_f,t_0,t_f+1e-6,structure.uniform_els,1,0,[],[]);
    dfdt = (f2-c)'/1e-6;
    Jc = [dfdt;Jc];
    
else
   
    c  = [];
    Jc = [];
    
end

[ceq,Jceq] = eval_constraints(f,structure,x,x_0,x_f,t_0,t_f,structure.uniform_els,1,dfx,dfu);

Jceq = Jceq';

% min time formulation
f2 = eval_constraints(f,structure,x,x_0,x_f,t_0,t_f+1e-6,structure.uniform_els,0,[],[]);
dfdt = (f2-ceq)'/1e-6;
Jceq = [dfdt;Jceq];

end