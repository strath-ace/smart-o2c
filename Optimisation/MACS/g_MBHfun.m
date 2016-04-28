function g=g_MBHfun(x,func,lambda,z,arg)

f=func(x,arg{:});

g=lambda*((f-z).^2)';

return