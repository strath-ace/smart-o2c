% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function x = lhsgen(n,p)

% generates latin hypesquare sample of cardinality n in p dimensions with ranges [0,1]
    for i=1:p
        x(:,i) = randperm(n);
    end
    x = (x - rand(n,p))/n;
    
end