% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [g,f,c,y]=ZDtestfront(test,z)
%
%   [f,c,y]=ZDtestfront(test,z)
%
%   INPUT
%               test : test function name
%                z     : test set
%
%   OUTPUT
%               g : multiobjective vector
%

% (c) Massimiliano Vasile

if isempty(z)
%x(:,1)=cos([0:0.005:pi/2]');
%x(:,2)=sin([0:0.005:pi/2]');
x(:,1)=0:0.01:1;
x(:,2:30)=zeros(length(x(:,1)),29);
%x(:,1)=1:0.01:5;
else
    x(:,:)=z(:,:);
end
lx=length(x(:,1));
c=[];
g=[];
length(test)
if length(test)>4
    for i=1:lx,
        for j=1:lx
            y(i,:)=[x(i,1) x(j,2) x(i,3:22)*0];
            [f(j+lx*(i-1),:),c(i,:)]=ZDtestfun(y(i,:),test);
        end
    end;
dom=dominance(f,0);
j=0;
for i=1:lx*lx,
    if dom(i)==0
        j=j+1;
        g(j,:)=f(i,:);
    else
        g(j,:)=[NaN NaN];
    end
end
    
else
     
    for i=1:lx,
        y(i,:)=x(i,:); %[x(i,1) x(i,2:end)*0];
        [f(i,:),c(i,:)]=ZDtestfun(y(i,:),test);
       
    end;
dom=dominance(f,0);
j=0;
for i=1:lx,
    if dom(i)==0
        j=j+1;
        g(j,:)=f(i,:);
    else
        j=j+1;
        g(j,:)=[NaN NaN];
        
    end
end
   
end


return