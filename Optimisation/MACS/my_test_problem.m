function [f]=my_test_problem(x,bogus)

    f(1) = sin(x(1)+x(2));
    f(2) = cos(x(1)+x(2));
    f(3) = x(1)+x(2);

end