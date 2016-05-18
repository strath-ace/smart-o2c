function [c,ceq]=c_MBHfun(x,func,lambda,z,arg)

f=func(x,arg{:});
c(1)=0;
ceq(1)=0;
lf=length(f);
index=1:lf;
for i=1:lf-1
    c(i)=lambda(index(i+1))*(f(index(i))-z(index(i)))-lambda(index(i))*(f(index(i+1))-z(index(i+1)));
end
ceq=c;
return
