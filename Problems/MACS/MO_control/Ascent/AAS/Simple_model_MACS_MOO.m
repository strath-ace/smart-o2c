function [val,x_sol] = Simple_model_MACS_MOO(x_in,lb,ub,structure,x_0,x_f,fminconoptions)

%fminconoptions = optimset('Algorithm','sqp');

x_in = (x_in-structure.offset_optimisation_vars)./structure.scale_optimisation_vars;
x_0 = (x_0-structure.offset_primal_states)./structure.scale_primal_states;
x_f = (x_f-structure.offset_primal_states)./structure.scale_primal_states;
lb = structure.norm_vlb;
ub = structure.norm_vub;

t_0 = 0;
t_f = x_in(1);

jacflag = 0;

%[x_sol,~,exitflag,~,~] = fmincon(@(x) 1,x_in(2:end),[],[],[],[],lb(2:end),ub(2:end),@(x) dynamics(structure,x,x_0,x_f,t_0,t_f,jacflag),fminconoptions);
[x_sol,~,exitflag,~,~] = fmincon(@(x) 1,x_in(2:end),[],[],[],[],lb(2:end),ub(2:end),@(x) dynamics(structure,x,x_0,x_f,t_0,t_f,jacflag),fminconoptions);

x_sol = [t_f x_sol];

x_sol = x_sol.*structure.scale_optimisation_vars+structure.offset_optimisation_vars;

t_f = x_sol(1);

if exitflag==1
             
    val = eval_control_objectives(x_sol(2:end),structure,x_f,t_0,t_f);
    
    % very ugly manual denormalisation of objective values
    
    %val(1) = val(1)*structure.scale_primal_states(5)+structure.offset_primal_states(5);
    %val(2) = val(2)*structure.scale_tf;
    
else
        
    val = [0 650];
    
end

end