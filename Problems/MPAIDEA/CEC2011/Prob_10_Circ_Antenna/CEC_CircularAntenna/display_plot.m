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
function[]=display_plot(gbest,phi_desired,distance)

dim=length(gbest);
phi=linspace(0,360,5000);
yax(1)=array_factorcir(gbest,(pi/180)*phi(1),phi_desired,distance,dim);
maxi=yax(1);
for i=2:5000%This loop finds out the maximum gain 
    yax(i)=array_factorcir(gbest,(pi/180)*phi(i),phi_desired,distance,dim);
    if maxi<yax(i)
        maxi=yax(i);
    end;
end;
for i=1:5000%This loop normalizes the Y-axis and finds the normalized gain values in decibels 
    yax(i)=yax(i)/maxi;
    yax(i)=20*log10(yax(i));
end;
plot(phi,yax,'g')
xlabel('Azimuth angle(deg)');
ylabel('Gain(db)');