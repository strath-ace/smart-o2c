function [out] = TC_1_constraints(d, u, par)


if u(1) < 0
    out = abs(u + sum(d));
else
    out = 0;
end

out = sum(out);
return

