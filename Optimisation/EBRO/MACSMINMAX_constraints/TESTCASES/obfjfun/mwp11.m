function [f] = mwp11(d,u,par)

f = cos(sqrt(d(1)^2+u(1)^2))/(10+sqrt(d(1)^2+u(1)^2));
global nfevalglobal
nfevalglobal = nfevalglobal + 1;

end