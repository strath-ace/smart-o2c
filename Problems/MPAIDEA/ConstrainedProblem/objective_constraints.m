function [y, c, ceq] = objective_constraints(x)
y = 5 .* (x - 3)^2 + 2;
c = 2 - x;
ceq = [];
end