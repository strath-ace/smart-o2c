function output = isPositiveIntegerValuedNumeric(x)

if all(floor(x) == x) && all(x > 0) 
    output = 1;
else
    output = 0;
end

return

