function [f] = mv3(d,u,par)
% Multidimensional test function MV3

f = sum( (5-d).*(1+cos(u)) + (d-1).*(1+sin(u)) );

global nfevalglobal
nfevalglobal = nfevalglobal + 1;

end
