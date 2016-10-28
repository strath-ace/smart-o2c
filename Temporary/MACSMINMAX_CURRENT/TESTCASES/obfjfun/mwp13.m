function [f] = mwp13(d,u,par)

f = (d(1)-2)^2+(d(2)-1)^2+u(1)*(d(1)^2-d(2))+u(2)*(d(1)+d(2)-2);
global nfevalglobal
nfevalglobal = nfevalglobal + 1;

end