function [f] = mv9(d,u,par)
% Multidimensional test function MV9

f = sum((d-u).*cos(-5*u+3*d));
global nfevalglobal
nfevalglobal = nfevalglobal + 1;

end
