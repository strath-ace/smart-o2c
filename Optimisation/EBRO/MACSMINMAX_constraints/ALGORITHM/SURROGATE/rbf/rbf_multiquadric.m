function psi = rbf_multiquadric(r,s)

if isempty(s), s = 1; end

psi = (r.^2 + s^2).^0.5;

end
