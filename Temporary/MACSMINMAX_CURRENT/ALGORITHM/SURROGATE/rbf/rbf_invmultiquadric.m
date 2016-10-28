function psi = rbf_invmultiquadric(r,s)

if isempty(s), s = 1; end

psi = (r.^2 + s^2).^-0.5;

end
