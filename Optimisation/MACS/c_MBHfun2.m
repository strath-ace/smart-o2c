function [c,ceq]=c_MBHfun2(x,func,lambda,f0,arg)

f=func(x,arg{:});
c = (lambda).*(f-f0);
ceq = [];

return
