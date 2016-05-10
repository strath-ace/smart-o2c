function val = objectives_rise (x_in,t_0,x_f,structure)

t_f = 80;
x_sol = x_in;

g = @(x,u,t) [-1/2*(x(2)^2+x(4)^2)+1/x(1) 0];
weights = [1 0];

[x,u,x_b] = extract_solution(x_sol,structure,x_f);

val = eval_cost_functions(g,weights,x,u,x_b,[t_0 t_f],structure,0,[],[],[],[]);

end