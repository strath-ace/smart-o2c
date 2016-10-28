function [f] = mwp8(d,u,par)

f = (d(1)-5)^2-(u(1)-5)^2;
global nfevalglobal
nfevalglobal = nfevalglobal + 1;

end