<<<<<<< HEAD
% Downloaded from:
% http://web.mysites.ntu.edu.sg/epnsugan/PublicSite/Shared%20Documents/Forms/AllItems.aspx?RootFolder=%2fepnsugan%2fPublicSite%2fShared%20Documents%2fCEC2005&FolderCTID=&View=%7bDAF31868%2d97D8%2d4779%2dAE49%2d9CEC4DC3F310%7d
=======
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%-----------Copyright (C) 2018 University of Strathclyde and Authors-----------
%

>>>>>>> 56eeb4b328a0319b2f58a2e2413248a83fcc168b

% function function_plot
% clear;close all
global initial_flag

% Xmin=[-100,-100,-100,-100,-100,-100,0,-32,-5,-5,-0.5,-pi,-3,-100,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,2];
% Xmax=[100,100,100,100,100,100,600,32,5,5,0.5,pi,1,100,5,5,5,5,5,5,5,5,5,5,5];

for func_num=1:1
    if func_num==1 x=-100:5:100;y=x; %[-100,100]
    elseif func_num==2 x=-100:5:100; y=x;%[-10,10]
    elseif func_num==3 x=-100:5:100; y=x;%[-100,100]
    elseif func_num==4 x=-100:2:0; y=x+100;%[-100,100]
    elseif func_num==5 x=-200:20:200;y=x; %[-5,5]
    elseif func_num==6 x=78:0.05:82; y=-51.3:0.05:-47.3;%[-3,1]
    elseif func_num==7 x=-350:2:-250; y=-100:2:0;;%[-600,600]
    elseif func_num==8 x=-32:1:32; y=x;%[-32,32]
    elseif func_num==9 x=-5:0.1:5; y=x;%[-5,5]
    elseif func_num==10 x=-5:0.1:5; y=x;%[-5,5]
    elseif func_num==11 x=-0.5:0.01:0.5; y=x;%[-0.5,0.5]
    elseif func_num==12 x=-3:0.1:3; y=x;%[-pi,pi]
    elseif func_num==13 x=-2:0.02:-1; y=x;%[-3,1]
    elseif func_num==14 x=-90:0.2:-50; y=-40:0.2:0;%[-100,100]
    elseif func_num==15 x=-5:0.1:5; y=x;%[-5,5]
    elseif func_num==16 x=-5:0.1:5; y=x;%[-5,5]
    elseif func_num==17 x=-5:0.1:5; y=x;%[-5,5]
    elseif func_num==18 x=-5:0.1:5; y=x;%[-5,5]
    elseif func_num==19 x=-5:0.1:5; y=x;%[-5,5]
    elseif func_num==20 x=-5:0.1:5; y=x;%[-5,5]        
    elseif func_num==21 x=-5:0.1:5; y=x;%[-5,5]
    elseif func_num==22 x=-5:0.1:5; y=x;%[-5,5]     
    elseif func_num==23 x=-5:0.1:5; y=x;%[-5,5]
    elseif func_num==24 x=-5:0.1:5; y=x;%[-5,5]   
    elseif func_num==25 x=-5:0.1:5; y=x;%[-5,5]  
    end
    
initial_flag=0;
L=length(x);
f=[];

for i=1:L
for j=1:L
f(i,j)=benchmark_func([x(i),y(j)],func_num);
end
end

figure(func_num)
surfc(x,y,f);
end

