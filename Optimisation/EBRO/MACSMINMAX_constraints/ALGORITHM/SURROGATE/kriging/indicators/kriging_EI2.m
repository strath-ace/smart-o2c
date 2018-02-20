% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function indicator = kriging_EI2(y, mse, ymin)

% [y, mse] = surrogate.predictor(x, surrogate.model);

s = sqrt(abs(mse));
%s(s <= 1e-4) = 0;

% Expected improvement
I = ymin - y;
% if all(I <= 0) || any(s == 0)
%     E = zeros(size(I));
% else
    u = I./s;
    PHI = 0.5*erfc(-u/sqrt(2)); % PHI = normcdf(u);
    phi = exp(-0.5*u.^2) ./ sqrt(2*pi); % phi = normpdf(u);
    E = I.*PHI+s.*phi+I.*exp(-s);
    % E(I==0 || s==0) = 0; % Intends to make it fail-safe for degenerate cases, I think it's OK but not sure... maybe should be I<=0?
    
% end

indicator = E;
        

end