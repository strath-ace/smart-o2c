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