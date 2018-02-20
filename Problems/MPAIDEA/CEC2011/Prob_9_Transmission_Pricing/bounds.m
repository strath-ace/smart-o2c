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
EBEinputfile; 

% n=length(bus_spec(:,1));
Pg=(bus_spec(:,7))/100;
Pd=(bus_spec(:,5))/100;

g=find(Pg>0);
d=find(Pd>0); 

%%%%%%%%% define BT   %%%%%%%%%%%%%%%%%%%%%
BT=zeros(length(g),length(d));

BT(1,4)=5;BT(1,5)=10;BT(1,6)=5;
BT(2,3)=5;
BT(3,21)=2.5;
BT(4,21)=2.5;BT(4,16)=15;
BT(5,12)=2.5;BT(6,8)=2.5;

BT=BT/100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

GD_max=zeros(length(g),length(d));
for i=1:length(g)
    for j=1:length(d)
        GD_max(i,j)=min(Pg(g(i))-BT(i,j),Pd(d(j))-BT(i,j));
    end
end



