function [c,ceq,Jc,Jceq] = dynamics_rise (f,structure,x_in,x_0,x_f,t_0,dfx,dfu)

t_f = 80;
x = x_in;
c = [];
[ceq] = eval_constraints(f,structure,x,x_0,x_f,t_0,t_f,0,[],[]);

Jc = [];

[ceq,Jceq] = eval_constraints(f,structure,x,x_0,x_f,t_0,t_f,1,dfx,dfu);

Jceq = Jceq';

end