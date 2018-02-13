function [minmax] = evaluate_minmax(problem_minmax, algo_minmax, algo_outer, algo_inner)






problem_minmax.sign_inner = 1;  %  for minmax



[ dmin, fminmax, exitflag, output ] = algo_minmax.optimise(problem_minmax,algo_outer,algo_inner,algo_minmax.par_minmax); % algo_minmax.optimise



minmax.d = dmin;
minmax.u = output.u;
minmax.f = fminmax;

end