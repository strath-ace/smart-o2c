function [c,ceq,Jc,Jceq] = Simple_dynamics (f,structure,x_in,x_0,x_f,t_0,dfx,dfu)

t_f = x_in(1);
x = x_in(2:end);

c = [];
Jc = [];

[ceq,Jceq] = eval_constraints(f,structure,x,x_0,x_f,t_0,t_f,structure.uniform_els,1,dfx,dfu);

%ceq = ceq.*structure.scale_transcribed_vars(2:end)/structure.scale_transcribed_vars(1);

Jceq = Jceq';

f2 = eval_constraints(f,structure,x,x_0,x_f,t_0,t_f+1e-6,structure.uniform_els,0,[],[]);
dfdt = (f2-ceq)'/1e-6;
Jceq = [dfdt;Jceq];

end