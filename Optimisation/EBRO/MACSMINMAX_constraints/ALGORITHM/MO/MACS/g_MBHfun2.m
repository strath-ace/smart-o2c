function g=g_MBHfun2(x,func,lambda,z,zstar,arg)

f=func(x,arg{:});

%g=lambda*((f-z).^2)';
%f2 = (f-z)./((zstar-z));

% in case zstar==z we have a Nan... to avoid it, i don't normalise in that
% case (also because in the single objective case zstar always equals z)

f2 = (f-z)./(((zstar-z).*((zstar-z)~=0))+((zstar-z)==0));

% trick to keep the sign :D
[~,index]= max(lambda.*abs(f2));
g = lambda(index)*f2(index);

%g = f(index);

return