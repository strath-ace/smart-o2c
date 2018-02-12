function [minmin] = evaluate_minmin(problem_minmax, algo_minmax, algo_outer, algo_inner)





problem_minmax.sign_inner = -1;  %  for minmin


[ dmin, fminmin, exitflag, output ] = algo_minmax.optimise(problem_minmax,algo_outer,algo_inner,algo_minmax.par_minmax); % algo_minmax.optimise


minmin.d = dmin;
minmin.u = output.u;
minmin.f = fminmin;

end