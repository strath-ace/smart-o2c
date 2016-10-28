function [masked] = mask_objfun_macsminmax_outer(d,par_objfun)

[masked, ~, ~] = u_validation( par_objfun.problem_fix_d, d, par_objfun.u_record, par_objfun.local_search, par_objfun.objectives);

return
