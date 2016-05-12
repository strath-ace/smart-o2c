function [c,ceq,Jc,Jceq] = dynamics2(f,structure,x_in,x_0,x_f,t_f,dfx,dfu)

c = [];
Jc = [];

t_0 = 0;
%t_f = x_in(1);
x_guess = x_in;%(2:end);

[ceq,Jceq] = eval_constraints(f,structure,x_guess,x_0,x_f,t_0,t_f,1,dfx,dfu);

Jceq = Jceq';
%f2 = eval_constraints(f,structure,x_guess,x_0,x_f,t_0,t_f+1e-6,0,[],[]);
%dfdt = (f2-ceq)'/1e-6;
%Jceq = [dfdt;Jceq];

end