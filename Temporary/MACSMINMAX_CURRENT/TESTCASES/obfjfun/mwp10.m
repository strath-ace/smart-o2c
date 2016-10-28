function [f] = mwp10(d,u,par)

f = sin(d(1)-u(1))/sqrt(d(1)^2+u(1)^2);
global nfevalglobal
nfevalglobal = nfevalglobal + 1;

end