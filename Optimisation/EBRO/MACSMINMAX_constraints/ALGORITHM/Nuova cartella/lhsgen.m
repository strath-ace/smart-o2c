function x = lhsgen(n,p)

% generates latin hypesquare sample of cardinality n in p dimensions with ranges [0,1]
    for i=1:p
        x(:,i) = randperm(n);
    end
    x = (x - rand(n,p))/n;
    
end