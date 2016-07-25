% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%-----------Copyright (C) 2016 University of Strathclyde-------------
%


function g=hp_g_fun(f,lambda,z,zstar)

% High performance g_fun. Identical to original one but is vectorized,
% sparing a whole for loop. Works on a single subproblem (lambda) for now
%
% INPUT
%       f       :       vector of objective function values
%       lambda  :       weigths vector of current problem
%       z       :       minimas of each objective function
%
% OUTPUT
%       g       :       vector of g_fun values

n_a = size(f,1);                                                            % number of agents

n_l = size(lambda,1);                                                       % number of subproblems

%g = max(repmat(lambda,n_a,1).*abs(f-repmat(z,n_a,1)),[],2)';                % efficient evaluation of g_function

g = zeros(n_a,n_l);
f2 = (f-repmat(z,n_a,1))./repmat(zstar-z,n_a,1);

for i = 1:n_l

    g(:,i) = max(repmat(lambda(i,:),n_a,1).*abs(f2),[],2)';

end
   

%d1 = zeros(n_a,1);
%d2 = zeros(n_a,1);

%for i=1:n_a

%    d1(i) = abs((f(i,:)-z)*lambda');
%    d2(i) = ((f(i,:)-(z+d1(i)*lambda))*(f(i,:)-(z+d1(i)*lambda))')^0.5;

%end

%g = d1+5*d2;

return