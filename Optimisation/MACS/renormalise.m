function [dd,energy,ener2,mins,maxs]=renormalise(memories,mins,maxs,lx,mfit)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%
% Renormalises the archive with the new max and min

f = memories(:,lx+1:lx+mfit);

delta = maxs-mins;
delta = delta.*(delta~=zeros(size(delta)))+ones(size(delta)).*(delta==zeros(size(delta))); % if there's only one value, max and min are the same, so this avoids a useless division by zero

f = f-repmat(mins,size(f,1),1);
f = f./repmat(delta,size(f,1),1);

xx = permute(f,[1,3,2]);                                                    % FAST computation of all relative square distances of candidates
xx = repmat(xx,1,size(f,1));
xxt = permute(xx,[2,1,3]);
dd = sum((xx-xxt).^2,3);

dd = dd.*(dd>0)+realmin.*(dd==0);                                           % Computation of pairwise energies, bottom right diagonal block
dd = 1./dd;

energy = sum(sum(triu(dd,1)));                                              % Energy

ener2 = zeros(size(memories,1),1);

for i = 1:size(memories,1)
    
    ener2(i) = energy-sum(dd(i,1:end~=i));
    
end


end