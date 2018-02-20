% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function s=lhsu(xmin,xmax,nsample,varargin)
% Latin Hypercube Sampling from uniform distribution
% Input:
%   xmin    : min of data (1,nvar)
%   xmax    : max of data (1,nvar)
%   nsample : no. of samples
% Output:
%   s       : random sample (nsample,nvar)
%   Budiman (2003)

nvar=length(xmin);
ran=rand(nsample,nvar);
s=zeros(nsample,nvar);

if ~isempty(varargin)
   
    ids = varargin{1};
    
else
    
    ids = 1:nvar;
    
end
for j=ids
   idx=randperm(nsample);
   P =(idx'-ran(:,j))/nsample;
   s(:,j) = xmin(j) + P.* (xmax(j)-xmin(j));
end