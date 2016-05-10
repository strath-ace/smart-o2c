function val = eval_control_objectives(x,structure,x_f,t_0,t_f)

[x_out,u,x_b] = extract_solution(x,structure,x_f);       
     
val = eval_cost_functions(structure.g,structure.weights,x_out,u,x_b,[t_0 t_f],structure,0,[],[])';

end