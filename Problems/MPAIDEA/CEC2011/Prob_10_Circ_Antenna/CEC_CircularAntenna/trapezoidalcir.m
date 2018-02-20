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
function[q]=trapezoidalcir(x2,upper,lower,N1,phi_desired,distance,dim)         
%This function performs integration by trapezoidal rule
h=(upper-lower)/N1;
 x1=lower;
 y=abs(realpow(array_factorcir(x2,lower,phi_desired,distance,dim),2)*sin(lower-pi/2));
%  y=abs(realpow(array_factorcir(x2,lower),2));
 for i=2:N1+1
 
  x1=x1+h;
  y(i)=abs(realpow(array_factorcir(x2,x1,phi_desired,distance,dim),2)*sin(x1-pi/2));
%   y(i)=abs(realpow(array_factorcir(x2,x1),2));
 end;
s=0;
 for i=1:N1+1
 
  if(i==1||i==N1+1)
	s=s+y(i);
  else
	s=s+2*y(i);
  end;
 end;
 q=(h/2)*s;

 
 




  




