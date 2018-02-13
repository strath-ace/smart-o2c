function [f] = mv8(d,u,par)
% Multidimensional test function MV8

f = sum((2*pi - u).*cos(u-d));
global nfevalglobal
nfevalglobal = nfevalglobal + 1;

end
