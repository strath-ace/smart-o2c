function [val,x_sol] = fastest_climb_MACS_fmincon2(x_in,lb,ub,f,structure,x_0,x_f,dfx,dfu,fminconoptions)

[x_sol,~,exitflag] = fmincon(@(x) 1,x_in,[],[],[],[],lb,ub,@(x) dynamics2(f,structure,x,x_0,x_f,dfx,dfu),fminconoptions);

t_f = x_sol(1);
t_0 = 0;

g = @(x,u,t) [t_f 0];
weights = [1 0];

[x_out,u,x_b] = extract_solution(x_sol(2:end),structure,x_f);

if exitflag==1
    
    val = eval_cost_functions(g,weights,x_out,u,x_b,[t_0 t_f],structure,0,[],[]);
    
else
    
    val = inf;
    
end

end

function [c,ceq,Jc,Jceq] = dynamics2(f,structure,x_in,x_0,x_f,dfx,dfu)

c = [];
Jc = [];

t_0 = 0;
t_f = x_in(1);
x_guess = x_in(2:end);

[ceq,Jceq] = eval_constraints(f,structure,x_guess,x_0,x_f,t_0,t_f,1,dfx,dfu);

Jceq = Jceq';
f2 = eval_constraints(f,structure,x_guess,x_0,x_f,t_0,t_f+1e-6,0,[],[]);
dfdt = (f2-ceq)'/1e-6;
Jceq = [dfdt;Jceq];

end