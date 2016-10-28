function indicator = kriging_EI(y, mse, ymin)

% [y, mse] = surrogate.predictor(x, surrogate.model);

s = sqrt(abs(mse));
%s(s <= 1e-4) = 0;

% Expected improvement
I = ymin - y;
% if all(I <= 0) || any(s == 0)
%     E = zeros(size(I));
% else
    u = I./s;
    PHI = normcdf(u);
    phi = normpdf(u);
    E = I.*PHI+s.*phi;
    % E(I==0 || s==0) = 0; % Intends to make it fail-safe for degenerate cases, I think it's OK but not sure... maybe should be I<=0?
    
% end

indicator = E;
        
end