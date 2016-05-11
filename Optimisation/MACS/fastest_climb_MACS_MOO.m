function [val,x_sol] = fastest_climb_MACS_MOO(x_in,lb,ub,f,structure,x_0,x_f,dfx,dfu,dft,fminconoptions)

[x_sol,~,exitflag,~,~] = fmincon(@(x) 1,x_in(2:end),[],[],[],[],lb(2:end),ub(2:end),@(x) dynamics2(f,structure,x,x_0,x_f,x_in(1),dfx,dfu),fminconoptions);

x_sol = [x_in(1) x_sol];

if exitflag==1
    
    t_f = x_sol(1);
    t_0 = 0;
    
    val = eval_control_objectives(x_sol(2:end),structure,x_f,t_0,t_f);

    %g = @(x,u,t) [t 0; -x(2) 0];    

    %weights = [1 0; 1 0];
    %[x_out,u,x_b] = extract_solution(x_sol(2:end),structure,x_f);       
     
    %val = eval_cost_functions(g,weights,x_out,u,x_b,[t_0 t_f],structure,0,[],[])';
        
else
        
    val = [350 5];

end

end