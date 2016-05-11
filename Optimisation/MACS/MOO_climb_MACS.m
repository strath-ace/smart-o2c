function [val,x_out] = MOO_climb_MACS(x_in,f,x_0,x_f,t_0,structure,tol_conv,maxits,dfx,dfu,lb,ub)

t_f = x_in(1);      % t_f is an optimisation variable!!!

x_guess = make_first_guess(f,x_0,t_0,x_in(1),x_in(2:end)',structure);

[x_sol,resids,t]=solve_constraints(f,x_guess,x_0,x_f,t_0,t_f,structure,tol_conv,maxits,dfx,dfu,lb,ub);    % solution of underdetermined system

g = @(x,u,t) [t_f 0; 0 u(1)]; %minimise indipendently climb time and integral of thrust over time, i.e. impulse ~ fuel consumption
weights = [1 0; 0 1];

num_obj = size(g(zeros(structure.num_eqs,1),zeros(structure.num_controls,1),0),1);

[x,u,x_b] = extract_solution(x_sol,structure,x_f);
x_out = [t_f u(:)'];

if resids(end)<tol_conv
    
    val = eval_cost_functions(g,weights,x,u,x_b,[t_0 t_f],structure,0,[],[]);
    
    
else

    val = inf(1,num_obj);

end