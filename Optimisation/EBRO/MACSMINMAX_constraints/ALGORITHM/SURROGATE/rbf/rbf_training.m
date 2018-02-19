% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function model = rbf_training(x, y, surrogate)

basis = surrogate.corrfun;
s = surrogate.s;

N = size(x,1);
R = zeros(N);

rbf_func = strcat('rbf_',lower(basis));
rbf_func = str2func(rbf_func);

for i = 1:N
    for j = 1:N
        R(i,j) = norm(x(i,:) - x(j,:));
    end
end
PSI = rbf_func(R,s);

w = PSI\y;

model.X = x;
model.Y = y;
model.PSI = PSI;
model.w = w;
model.basis = basis;
model.s = s;

end
