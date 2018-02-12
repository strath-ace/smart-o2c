function [masked, ceq] = mask_constraints_macsminmax_outer(d,par_objfun)   % par_objfun.u_record

%par_objfun.problem_fix_d.par_objfun.objfun = par_objfun.problem_fix_d.par_objfun.constraint;

%[masked, ~, ~] = u_validation( par_objfun.problem_fix_d, d, par_objfun.u_record, par_objfun.local_search, par_objfun.objectives);



d_true = par_objfun.problem_fix_d.par_objfun.lb_d' + d.*(par_objfun.problem_fix_d.par_objfun.ub_d' - par_objfun.problem_fix_d.par_objfun.lb_d');
obj = par_objfun.objectives;
map_u_info = par_objfun.problem_fix_d.par_objfun.map_u_info{obj};

% check the constraint in all the u in the archive

for i = 1:size(par_objfun.u_record{1},1)
    
    u_true = map_affine(par_objfun.u_record{1}(i,:), map_u_info);
    func = par_objfun.problem_fix_d.par_objfun.constraint{1}; 
    % func =  par_objfun.constraint{1};
    par_func = par_objfun.problem_par_objfun{obj};%  par_objfun.problem_fix_d.par_objfun;

    Constr(i) = par_objfun.problem_fix_d.par_objfun.sign*func(d_true,u_true,par_func);
    
end


masked = max(Constr);
ceq=[];
return
