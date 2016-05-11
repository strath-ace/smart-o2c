function [c,ceq,Jc,Jceq] = dynamics(structure,x_in,x_0,x_f,t_0,t_f,jacflag)

c = [];
Jc = [];

[ceq,Jceq] = eval_constraints(structure.f,structure,x_in,x_0,x_f,t_0,t_f,jacflag,structure.dfx,structure.dfu);

Jceq = Jceq';

end