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
%tersoff potential based on Si(B) model.
function f=tersoff(x)
R1=3.0;
R2=0.2;
A=3.2647e+3;
B=9.5373e+1;
lemda1=3.2394;
lemda2=1.3258;
lemda3=1.3258;
c=4.8381;
d=2.0417;
n1=22.956;
gama=0.33675;
h=0;
E=zeros(NP);
for i=1:NP
    for j=1:NP
         r(i,j)=sqrt(sum((x(i,:)-x(j,:)).^2)/3);
if r(i,j)<(R1-R2)
    fcr(i,j)=1;
elseif  r(i,j)>(R1+R2)
    fcr(i,j)=0;
else
    fcr(i,j)=0.5-0.5*sin(pi/2*(r(i,j)-R1)/R2);
end

VRr(i,j)=A*exp(-lemda1*r(i,j));
VAr(i,j)=B*exp(-lemda2*r(i,j));
jeta=zeros(NP,NP);
for k=1:NP
    if i==k || j==k        continue  
    end
     rd1=sqrt(sum((x(i,:)-x(k,:)).^2)/3);
      rd3=sqrt(sum((x(k,:)-x(j,:)).^2)/3);
       rd2=sqrt(sum((x(i,:)-x(j,:)).^2)/3);
       ctheta_ijk=(rd1^2+rd2^2-rd3^3)/(2*rd1*rd2);
       G_th_ijk =1+(c^2)/(d^2)-(c^2)/(d^2+(h-ctheta_ijk)^2);
       jeta(i,j)=jeta(i,j)+fcr(i,k)*G_th_ijk*exp(lemda3^3*(r(i,j)-r(i,k))^3);
end
if i==j
    continue
end
Bij=(1+(gama*jeta(i,j))^n1)^(-0.5/n1);
E(i)=E(i)+fcr(i,j)*(VRr(i,j)-Bij*VAr(i,j))/2;
    end    
end
f=sum(E);





    
