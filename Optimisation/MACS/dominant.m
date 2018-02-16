function [domA,domB]=dominant(A,B)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%
%  Evaluates the (simple) dominance of  elements of array A with respect to
%  array B and vice versa.
%  Elements of both arrays should be already "self-dominant" 
%  among each array.

%   A(i,:)>B(j,:)
%      simple dominance is computed, i.e.
%      i is dominated by j if all the not equal components of j are better than i

[nA,mA]=size(A);
[nB,mB]=size(B);
domA = zeros(nA,1);
domB = zeros(nB,1);
if mA~=mB

    error('Vectors must have the same number of columns')
    
end


for i=1:nB
    rB = repmat(B(i,:),nA,1);                                               % matrix tiling, to make all combinations
    eq = repmat(all(A==rB,2),1,mA);                                         % tiled equality matrix
    
    domA = domA + sum(reshape(min((A>=rB).*~eq,[],2),nA,1),2);              % dominance of element rB on A
    domB(i) = sum(reshape(min((rB>=A).*~eq,[],2),1,nA)+nB*nA*reshape(min((rB==A).*eq,[],2),1,nA),2);    % dominance of A on element rB
end
return