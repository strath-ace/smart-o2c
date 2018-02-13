function [positions] = var2opt(i,problem_minmax)


positions = sum(problem_minmax.dim_u_i(1:i-1))+1 : sum(problem_minmax.dim_u_i(1:i)); 


end