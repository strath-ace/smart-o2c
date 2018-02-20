% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function  [f, df] = regpoly0(S)
%REGPOLY0  Zero order polynomial regression function
%
% Call:    f = regpoly0(S)
%          [f, df] = regpoly0(S)
%
% S  : m*n matrix with design sites
% f  : ones(m,1)
% df : Jacobian at the first point (first row in S) 

% hbn@imm.dtu.dk  
% Last update  April 12, 2002

[m n] = size(S);
f = ones(m,1);
if  nargout > 1
  df = zeros(n,1);
end
