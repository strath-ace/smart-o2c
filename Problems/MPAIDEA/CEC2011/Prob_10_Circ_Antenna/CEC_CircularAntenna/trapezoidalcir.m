% Downloaded from:
% http://web.mysites.ntu.edu.sg/epnsugan/PublicSite/Shared%20Documents/Forms/AllItems.aspx?RootFolder=%2fepnsugan%2fPublicSite%2fShared%20Documents%2fCEC%202011%2d%20RWP&FolderCTID=&View=%7bDAF31868%2d97D8%2d4779%2dAE49%2d9CEC4DC3F310%7d


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

 
 




  




