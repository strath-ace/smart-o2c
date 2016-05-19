function [c,ceq] = ascent_constr(x)

% x(1) = beta_0
% x(2) = beta_f
% x(3) = T_f

a = 4e-3;
g = 1.6e-3;
H = 10;

b0 = x(1);
bf = x(2);
Tf = x(3);

c = (tan(b0)-tan(bf))/Tf;

ceq = [ a/c*(sec(b0)-sec(bf))-g*Tf;
        a/(2*c^2)*( (tan(b0)-tan(bf))*sec(b0) - (sec(b0)-sec(bf)).*tan(bf) - log( (tan(b0)+sec(b0))./(tan(bf)+sec(bf))) )-0.5*g*Tf.^2-H];

c = [];

end