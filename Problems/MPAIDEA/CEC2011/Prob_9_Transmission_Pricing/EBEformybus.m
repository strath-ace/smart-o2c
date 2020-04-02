% Downloaded from:
% http://web.mysites.ntu.edu.sg/epnsugan/PublicSite/Shared%20Documents/Forms/AllItems.aspx?RootFolder=%2fepnsugan%2fPublicSite%2fShared%20Documents%2fCEC%202011%2d%20RWP&FolderCTID=&View=%7bDAF31868%2d97D8%2d4779%2dAE49%2d9CEC4DC3F310%7d

% Formation of Ybus.
function [YIbus] = formybus(linedata,n)
busa=linedata(:,1);busb=linedata(:,2);
L=length(busa);
Z=linedata(:,3);ZI=imag(Z);
YI=ones(L,1)./ZI;
YIbus=zeros(n,n);

    for k=1:L
    
    n=busa(k);m=busb(k);
    YIbus(n,n)=YIbus(n,n)+YI(k);
    YIbus(n,m)=-YI(k)+YIbus(n,m);
    YIbus(m,n)=-YI(k)+YIbus(m,n);
    YIbus(m,m)=YIbus(m,m)+YI(k);    
    end

return;
