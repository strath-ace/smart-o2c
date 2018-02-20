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

 
 




  




