% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function S = lhsamp(m, n)
%LHSAMP  Latin hypercube distributed random numbers
%
% Call:    S = lhsamp
%          S = lhsamp(m)
%          S = lhsamp(m, n)
%
% m : number of sample points to generate, if unspecified m = 1
% n : number of dimensions, if unspecified n = m
%
% S : the generated n dimensional m sample points chosen from
%     uniform distributions on m subdivions of the interval (0.0, 1.0)

% hbn@imm.dtu.dk  
% Last update April 12, 2002

if nargin < 1, m = 1; end
if nargin < 2, n = m; end

S = zeros(m,n);
for i = 1 : n
  S(:, i) = (rand(1, m) + (randperm(m) - 1))' / m;
end