function [c,ceq] = ascent_specific_smooth_scal_constraints(in,lambda,z,zstar,options)

t_0 = 0;
t_f = in(2);

[c,ceq]=smooth_cheb_constr(in,lambda,z,zstar,options,t_0,t_f,options.oc.x_0,options.oc.x_f);

end