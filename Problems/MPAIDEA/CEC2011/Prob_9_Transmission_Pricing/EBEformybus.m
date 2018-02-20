<<<<<<< HEAD
% Downloaded from:
% http://web.mysites.ntu.edu.sg/epnsugan/PublicSite/Shared%20Documents/Forms/AllItems.aspx?RootFolder=%2fepnsugan%2fPublicSite%2fShared%20Documents%2fCEC%202011%2d%20RWP&FolderCTID=&View=%7bDAF31868%2d97D8%2d4779%2dAE49%2d9CEC4DC3F310%7d

=======
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
>>>>>>> 56eeb4b328a0319b2f58a2e2413248a83fcc168b
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
