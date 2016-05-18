function [c,ceq,Jc,Jceq] = dynamics (f,structure,x_in,x_0,x_f,t_0,dfx,dfu)

% min time formulation
t_f = x_in(1);
x = x_in(2:end);

% max vel formulation
%t_f = 200;
%x = x_in;


c = [];
%[ceq] = eval_constraints(f,structure,x,x_0,x_f,t_0,t_f,0,[],[]);

Jc = [];

[ceq,Jceq] = eval_constraints(f,structure,x,x_0,x_f,t_0,t_f,1,dfx,dfu);

Jceq = Jceq';

% min time formulation
f2 = eval_constraints(f,structure,x,x_0,x_f,t_0,t_f+1e-6,0,[],[]);
dfdt = (f2-ceq)'/1e-6;
Jceq = [dfdt;Jceq];

end