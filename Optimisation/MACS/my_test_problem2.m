function [f]=my_test_problem2(x,bogus)

    f(1) = (x(1)+x(2)).^2;
    %f(2) = sin(5*pi*(x(1)-x(2)));
    f(2) = (x(1)-x(2)-1).^2;
    f(3) = (x(1)+2*x(2)-2).^2;

end