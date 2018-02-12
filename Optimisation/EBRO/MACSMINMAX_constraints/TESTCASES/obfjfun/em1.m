function [f] = em1(d,u,par)
% Multidimensional test function MV2

f = sum( (u-3*d).*sin(u)+(d-2).^2 );
global nfevalglobal
nfevalglobal = nfevalglobal + 1;

end
