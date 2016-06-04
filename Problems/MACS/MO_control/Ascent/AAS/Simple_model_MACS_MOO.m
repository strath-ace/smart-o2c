function [val,x_sol] = Simple_model_MACS_MOO(x_in,lb,ub,structure,x_0,x_f,fminconoptions)

%fminconoptions = optimset('Algorithm','sqp');

t_0 = 0;
t_f = x_in(1);

jacflag = 1;

[x_sol,~,exitflag,~,~] = fmincon(@(x) 1,x_in(2:end),[],[],[],[],lb(2:end),ub(2:end),@(x) dynamics(structure,x,x_0,x_f,t_0,t_f,jacflag),fminconoptions);

x_sol = [t_f x_sol];

if exitflag==1
             
    val = eval_control_objectives(x_sol(2:end),structure,x_f,t_0,t_f);
    
else
        
    val = [0 0];
    
end

end