function [f,df]=smooth_cheb_fun(x)

f = x(1);
df = 0*x;
df(1) = 1;

return