% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [masked, ceq] = mask_constraints_macsminmax_outer(d,par_objfun)   



d_true = par_objfun.problem_fix_d.par_objfun.lb_d' + d.*(par_objfun.problem_fix_d.par_objfun.ub_d' - par_objfun.problem_fix_d.par_objfun.lb_d');
obj = par_objfun.objectives;
map_u_info = par_objfun.problem_fix_d.par_objfun.map_u_info{obj};

% check the constraint in all the archive Au = Af U Ac

for i = 1:size(par_objfun.u_record{1},1)
    
    u_true = map_affine(par_objfun.u_record{1}(i,:), map_u_info);
    func = par_objfun.problem_fix_d.par_objfun.constraint{1}; 
    par_func = par_objfun.problem_par_objfun{obj};

    Constr(i) = par_objfun.problem_fix_d.par_objfun.sign*func(d_true,u_true,par_func);
    
end


masked = max(Constr);
ceq=[];
return
