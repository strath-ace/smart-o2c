<<<<<<< HEAD
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
=======
>>>>>>> 5b7361d93c9119cf1d2e9e6c885bed93f924d71b
function[y]=array_factorcir(x1,phi,phi_desired,distance,dim)
%As the function name signifies here the array factor is calculated 
 
pi=3.141592654;
  y=0;
  y1=0;
%   for i1=1:dim/2
%       delphi=2*pi*(i1-1)/dim;
%       shi=cos(phi-delphi)-cos(phi_desired*(pi/180)-delphi);
%       shi=shi*dim*distance;
%       y=y+x1(i1)*cos(shi+x1(dim/2+i1)*(pi/180));
%   end;
%   y=y*2;
%   y=abs(y);
 

  for i1=1:dim/2
      delphi=2*pi*(i1-1)/dim;
      shi=cos(phi-delphi)-cos(phi_desired*(pi/180)-delphi);
      shi=shi*dim*distance;
      y=y+x1(i1)*cos(shi+x1(dim/2+i1)*(pi/180));
  end;
   for i1=dim/2+1:dim
      delphi=2*pi*(i1-1)/dim;
      shi=cos(phi-delphi)-cos(phi_desired*(pi/180)-delphi);
      shi=shi*dim*distance;
      y=y+x1(i1-dim/2)*cos(shi-x1(i1)*(pi/180));
  end;
  for i1=1:dim/2
      delphi=2*pi*(i1-1)/dim;
      shi=cos(phi-delphi)-cos(phi_desired*(pi/180)-delphi);
      shi=shi*dim*distance;
      y1=y1+x1(i1)*sin(shi+x1(dim/2+i1)*(pi/180));
  end;
  for i1=dim/2+1:dim
      delphi=2*pi*(i1-1)/dim;
      shi=cos(phi-delphi)-cos(phi_desired*(pi/180)-delphi);
      shi=shi*dim*distance;
      y1=y1+x1(i1-dim/2)*sin(shi-x1(i1)*(pi/180));
  end;
  y=y*y+y1*y1;
  y=sqrt(y);
  