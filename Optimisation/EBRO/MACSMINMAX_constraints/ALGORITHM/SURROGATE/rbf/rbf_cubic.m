function [psi,dpsi] = rbf_cubic(r,s)

psi = r.^3;
dpsi = 3*r.^2;

end
