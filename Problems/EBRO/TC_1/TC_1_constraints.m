function [out] = TC_1_constraints(d, u, par)


if sum(u) >0
    out = u + sum(d);
else
    out = -1 + u - sum(d);
end

out = sum(out);
return

