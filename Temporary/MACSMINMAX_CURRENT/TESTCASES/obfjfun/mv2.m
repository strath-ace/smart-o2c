function [f] = mv2(d,u,par)
% Multidimensional test function MV2

f = sum( (d-u).^2 );
global nfevalglobal
nfevalglobal = nfevalglobal + 1;

end
