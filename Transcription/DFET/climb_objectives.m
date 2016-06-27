function [val,grad] = climb_objectives (x_in,t_0,x_f,structure)

% min time formulation
t_f = x_in(1);
x_sol = x_in(2:end);
g = @(x,u,t) [t 0];

% max vel formulation
%t_f = 250;
%x_sol = x_in;
%g = @(x,u,t) [-x(2) 0];

weights = [1 0];

%[x,u,x_b] = extract_solution(x_sol,structure,x_f);

[val,grad] = eval_cost_functions2(g,weights,x_sol,x_f,[t_0 t_f],structure.uniform_els,structure,1,[],[]);

val2 = eval_cost_functions2(g,weights,x_sol,x_f,[t_0 t_f+1e-6],structure.uniform_els,structure,0,[],[]);
dgdt = (val2-val)/1e-6;
grad = [dgdt grad];


end