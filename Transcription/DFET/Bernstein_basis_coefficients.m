function coeffs = Bernstein_basis_coefficients(k,n)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%
% Generates the coefficients of the k-th Bernstein Polynomial of degree n.
% Coefficients are NOT the standard Bernstein coefficients, but are
% converted in the canonical form tau^k (where -1<=tau<=1 !!!)

P = pascal(n+1,1);
p = P(end,:);
C1 = p'.*P;
C2 = diag((1/2).^(n:-1:0));
C3 = fliplr(flipud(abs(pascal(n+1,1))));
M = flipud(C1*C2*C3);

coeffs = M(k,:);

end