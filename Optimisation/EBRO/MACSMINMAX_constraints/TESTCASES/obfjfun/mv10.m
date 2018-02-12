function [f] = mv10(d,u,par)
% Multidimensional test function MV10

f = sum( (d+u).*cos(-u*(5*norm(d)+5)+3*d) );
global nfevalglobal
nfevalglobal = nfevalglobal + 1;

end
