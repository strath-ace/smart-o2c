% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [domA,domB,strong_dom]=dominant(A,B)
%
%   dominant
%
%  [domA,domB]=dominant(A,B);
%
%  Evaluates the (simple) dominance of  elements of array A with respect to
%  array B and vice versa.
%  Elements of both arrays should be already "self-dominant" 
%  among each array.
%
%   A(i,:)>B(j,:)
%      simple dominance is computed, i.e.
%      i is dominated by j if all the not equal components of j are better than i
%
%   (c)  Federico Zuiani  2010
%

strong_dom=0;
[nA,mA]=size(A);
[nB,mB]=size(B);
if mA~=mB
    error('Vectors must have the same number of columns')
end

domA=zeros(nA,1);
domB=zeros(nB,1);

for i=1:nA
    B_dom=0;
    for j=1:nB
        if A(i,:)==B(j,:)
            domB(j)=domB(j)+nA*nB;
        else
            domA(i)=domA(i)+min(A(i,:)>=B(j,:));
            domB(j)=domB(j)+min(B(j,:)>=A(i,:));
            B_dom=max([B_dom all(B(j,:)<A(i,:))]);
        end
    end
    strong_dom=strong_dom+B_dom/nA;
end

return