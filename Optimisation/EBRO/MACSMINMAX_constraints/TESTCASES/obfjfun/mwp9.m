function [f] = mwp9(d,u,par)

f = min([3-0.2*d(1)+0.3*u(1), 3+0.2*d(1)-0.1*u(1)]);
global nfevalglobal
nfevalglobal = nfevalglobal + 1;

end