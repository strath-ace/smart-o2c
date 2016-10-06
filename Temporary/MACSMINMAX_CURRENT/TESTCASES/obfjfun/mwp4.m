function [f] = mwp4(d,u,par)

f = -sum((u-1).^2)+sum((d-1).^2)+u(3)*(d(2)-1)+u(1)*(d(1)-1)+u(2)*d(1)*d(2);
global nfevalglobal
nfevalglobal = nfevalglobal + 1;

end