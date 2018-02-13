function [f] = mwp5(d,u,par)

f = -(d(1)-1)*u(1)-(d(2)-2)*u(2)-(d(3)-1)*u(3)+2*d(1)^2+3*d(2)^2+d(3)^2-u(1)^2-u(2)^2-u(3)^2;
global nfevalglobal
nfevalglobal = nfevalglobal + 1;

end