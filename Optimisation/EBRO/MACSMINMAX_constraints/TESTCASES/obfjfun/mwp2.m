function [f] = mwp2(d,u,par)

f = 4*(d(1)-2)^2-2*u(1)^2+d(1)^2*u(1)-u(2)^2+2*d(2)^2*u(2);
global nfevalglobal
nfevalglobal = nfevalglobal + 1;

end