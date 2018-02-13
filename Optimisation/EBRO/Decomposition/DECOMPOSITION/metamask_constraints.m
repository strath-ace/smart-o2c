function [masked, ceq] = metamask_constraints(u,par_objfun)
% unscales d, u, evaluates the corresponding objfun and multiplies by sign

% if par_objfun.sign == -1
%     par_objfun.d = u(1:par_objfun.dim_d);
%     u = u(1+par_objfun.dim_d:end);
%     par_objfun.objective = 1;
% end


    d_true = par_objfun.d;
%     d_true = par_objfun.lb_d' + par_objfun.d.*(par_objfun.ub_d' - par_objfun.lb_d');
    obj = par_objfun.objective;
    map_u_info = par_objfun.map_u_info{obj};

    u_true = map_affine(u,map_u_info);


    func = par_objfun.constraint{obj};
    par_func = par_objfun.problem_par_objfun{obj};

    masked = par_objfun.sign*func(d_true,u_true,par_func);


ceq = [];


return