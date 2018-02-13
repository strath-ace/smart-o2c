function [f] = mwp6(d,u,par)

f = u(1)*(d(1)^2-d(2)+d(3)-d(4)+2) + u(2)*(-d(1)+2*d(2)^2-d(3)^2 ...
            +2*d(4)+1)+u(3)*(2*d(1)-d(2)+2*d(3)-d(4)^2+5)+5*d(1)^2+4*d(2)^2+3*d(3)^2+2*d(4)^2-...
            sum(u.^2);
global nfevalglobal
nfevalglobal = nfevalglobal + 1;

end