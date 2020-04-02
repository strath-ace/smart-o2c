% Downloaded from:
% http://web.mysites.ntu.edu.sg/epnsugan/PublicSite/Shared%20Documents/Forms/AllItems.aspx?RootFolder=%2fepnsugan%2fPublicSite%2fShared%20Documents%2fCEC%202011%2d%20RWP&FolderCTID=&View=%7bDAF31868%2d97D8%2d4779%2dAE49%2d9CEC4DC3F310%7d


tol=1.0e-08;
tspan=[0 0.78];
yo =[ 0.09 0.09 ]';
u=2;
options = odeset('RelTol',tol);
[T,Y] = ode45(@(t,y) rigid(t,y,u),tspan,yo,options);
j=sum(sum(Y.^2,2)+(0.1)*u*u)
