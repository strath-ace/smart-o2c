% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [f,c]=cec2009test(x,mof)

c=0;

switch mof
    
    case 'UP1'
        
        n = length(x);
        j1 = 3:2:n;
        j2 = 2:2:n;
        f(1) = x(1) + 2/length(j1) * sum((x(j1)-sin(6*pi*x(1)+j1*pi/n)).^2);
        f(2) = 1-x(1)^0.5 + 2/length(j2) * sum((x(j2)-sin(6*pi*x(1)+j2*pi/n)).^2);

%         [dim, num]  = size(x);
%         tmp         = zeros(dim,num);
%         tmp(2:dim,:)= (x(2:dim,:) - sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*repmat((2:dim)',[1,num]))).^2;
%         tmp1        = sum(tmp(3:2:dim,:));  % odd index
%         tmp2        = sum(tmp(2:2:dim,:));  % even index
%         f(1,:)      = x(1,:)             + 2.0*tmp1/size(3:2:dim,2);
%         f(2,:)      = 1.0 - sqrt(x(1,:)) + 2.0*tmp2/size(2:2:dim,2);

end




end
