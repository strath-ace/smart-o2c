function [val,x_sol] = fastest_climb_MACS_MOO(x_in,lb,ub,f,structure,x_0,x_f,dfx,dfu,dft,fminconoptions)

x_init = x_0;
x_init(5) = x_in(2);

[x_sol,~,exitflag,~,~] = fmincon(@(x) 1,x_in(3:end),[],[],[],[],lb(3:end),ub(3:end),@(x) dynamics2(f,structure,x,x_init,x_f,x_in(1),dfx,dfu),fminconoptions);

x_sol = [x_in(1:2) x_sol];

if exitflag==1
    
    t_f = x_sol(1);
    t_0 = 0;
    
    g = @(x,u,t) [t 0; -x(2) 0];    

    weights = [1 0; 1 0];
    [x_out,u,x_b] = extract_solution(x_sol(3:end),structure,x_f);    
     
    val = eval_cost_functions(g,weights,x_out,u,x_b,[t_0 t_f],structure,0,[],[])';
        
else
        
    val = [350 5];

end

end

function [c,ceq,Jc,Jceq] = dynamics2(f,structure,x_in,x_0,x_f,t_f,dfx,dfu)

c = [];
Jc = [];

t_0 = 0;
x_guess = x_in;

[ceq,Jceq] = eval_constraints(f,structure,x_guess,x_0,x_f,t_0,t_f,1,dfx,dfu);

Jceq = Jceq';

end