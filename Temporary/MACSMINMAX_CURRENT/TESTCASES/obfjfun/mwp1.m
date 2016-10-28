function [f] = mwp1(d,u,par)

f = 5*(d(1)^2+d(2)^2)-(u(1)^2+u(2)^2)+d(1)*(-u(1)+u(2)+5)+ d(2)*(u(1)-u(2) + 3);
global nfevalglobal
nfevalglobal = nfevalglobal + 1;

end