function [masked] = meta_mask_objfun_min_u(u_to_opt,par_objfun)
% unscales d, u, evaluates the corresponding objfun and multiplies by sign

% d_true = par_objfun.lb_d' + par_objfun.d.*(par_objfun.ub_d' - par_objfun.lb_d');
% DA VERIFICARE 


d_true = par_objfun.d;                                 % for the moment I don't change d

obj = par_objfun.objective;
map_u_info = par_objfun.map_u_info{obj};

% u_aux = par_objfun.umax;

% u_aux(par_objfun.vars_to_opt)=u_to_opt;
u_true = map_affine(u_to_opt,map_u_info);

%map_u_info = par_objfun.map_u_info{obj};

% u_aux = par_objfun.umax;
% 
% % u_aux(par_objfun.vars_to_opt)=u_to_opt;
% % u_true = map_affine(u_aux,map_u_info);
% 
% %u_true = map_affine(u_to_opt,map_u_info);
% u_aux(par_objfun.vars_to_opt)=u_to_opt;
% u_true = u_aux;


func = par_objfun.objfun{obj};
par_func = par_objfun.problem_par_objfun{1};

masked = -par_objfun.flag*func(d_true,u_true,par_func);    % -par_objfun.sign*   VERIFICARE 


end