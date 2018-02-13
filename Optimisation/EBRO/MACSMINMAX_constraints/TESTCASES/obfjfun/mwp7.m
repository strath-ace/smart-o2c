function [f] = mwp7(d,u,par)

f = 2*d(1)*d(5)+3*d(4)*d(2)+d(5)*d(3)+5*d(4)^2+5*d(5)^2-d(4)*(u(4)-u(5)-5)+d(5)*(u(4)-u(5)+3)+sum(u(1:3).*(d(1:3).^2-1))-sum(u.^2);
global nfevalglobal
nfevalglobal = nfevalglobal + 1;

end