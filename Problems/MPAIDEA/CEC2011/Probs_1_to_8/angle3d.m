% Downloaded from:
% http://web.mysites.ntu.edu.sg/epnsugan/PublicSite/Shared%20Documents/Forms/AllItems.aspx?RootFolder=%2fepnsugan%2fPublicSite%2fShared%20Documents%2fCEC%202011%2d%20RWP&FolderCTID=&View=%7bDAF31868%2d97D8%2d4779%2dAE49%2d9CEC4DC3F310%7d

function th = angle3d(x,j,i,k)
%
%
a =x(i,:)-x(j,:);
b =x(k,:)-x(j,:);
A =sqrt(a(1)^2+a(2)^2+a(3)^2);
B =sqrt(b(1)^2+b(2)^2+b(3)^2);
c =dot(a,b);
if A*B == 0
    th=pi;
else
    th =acos(c/(A*B));
end
