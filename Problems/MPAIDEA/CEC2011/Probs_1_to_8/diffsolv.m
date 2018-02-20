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
function dy = diffsolv(t,x,u)
load c_bifunc_data;% c(i,j) is saved here.
ml=[1 u u.^2 u.^3];
mlt=repmat(ml,10,1);
k=sum(c.*mlt,2);
dy = zeros(7,1);    % a column vector
dy(1) = -k(1)*x(1);
dy(2) = k(1)*x(1)-(k(2)+k(3))*x(2)+k(4)*x(5);
dy(3) = k(2)*x(2);
dy(4) = -k(6)*x(4)+k(5)*x(5);
dy(5) = k(3)*x(2)+k(6)*x(4)-(k(4)+k(5)+k(8)+k(9))*x(5)+k(7)*x(6)+k(10)*x(7);
dy(6) = k(8)*x(5)-k(7)*x(6);
dy(7) = k(9)*x(5)-k(10)*x(7);
