function [f] = mwp3(d,u,par)

f = d(1)^4*u(2)+2*d(1)^3*u(1)-d(2)^2*u(2)*(u(2)-3)-2*d(2)*(u(1)-3)^2;
global nfevalglobal
nfevalglobal = nfevalglobal + 1;

end