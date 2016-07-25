% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%-----------Copyright (C) 2016 University of Strathclyde-------------
%


% Formation of Ybus.
function [YIbus] = formybus(linedata,n)
busa=linedata(:,1);busb=linedata(:,2);
L=length(busa);
Z=linedata(:,3);ZI=imag(Z);
YI=ones(L,1)./ZI;
YIbus=zeros(n,n);

    for k=1:L
    
    n=busa(k);m=busb(k);
    YIbus(n,n)=YIbus(n,n)+YI(k);
    YIbus(n,m)=-YI(k)+YIbus(n,m);
    YIbus(m,n)=-YI(k)+YIbus(m,n);
    YIbus(m,m)=YIbus(m,m)+YI(k);    
    end

return;