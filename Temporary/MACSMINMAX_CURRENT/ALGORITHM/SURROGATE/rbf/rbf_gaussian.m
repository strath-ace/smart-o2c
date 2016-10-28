function psi = rbf_gaussian(r,s)

if isempty(s), s = 1; end

psi = exp(-r.^2./(2*s^2));

end
