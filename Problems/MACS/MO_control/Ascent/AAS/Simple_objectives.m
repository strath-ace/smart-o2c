function [val] = Simple_objectives (x_in,t_0,x_f,structure)

% min fuel formulation
t_f = x_in(1);
x_sol = x_in(2:end);
g = @(x,u,t) [-x(5) 0];

% max vel formulation
% t_f = x_in(1);
% x_sol = x_in(2:end);
% g = @(x,u,t) [-x(3) 0];

weights = [1 0];

[x,u,x_b] = extract_solution(x_sol,structure,x_f);

val = eval_cost_functions(g,weights,x,u,x_b,[t_0 t_f],structure,0,[],[]);

end