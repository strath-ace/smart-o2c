function [f] = mv1(d,u,par)
% Multidimensional test function MV1

f = sum( d.*u.^2 );

global nfevalglobal
nfevalglobal = nfevalglobal + 1;

end
