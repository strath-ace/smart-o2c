% Downloaded from:
% http://web.mysites.ntu.edu.sg/epnsugan/PublicSite/Shared%20Documents/Forms/AllItems.aspx?RootFolder=%2fepnsugan%2fPublicSite%2fShared%20Documents%2fCEC%202011%2d%20RWP&FolderCTID=&View=%7bDAF31868%2d97D8%2d4779%2dAE49%2d9CEC4DC3F310%7d


function dy = intgrl(t,x,u)
dy = zeros(2,1);    % a column vector
dy(1) = -(2+u)*(x(1)+0.25)+(x(2)+0.5)*exp(25*x(1)/(x(1)+2));
dy(2) = 0.5-x(2)-(x(2)+0.5)*exp(25*x(1)/(x(1)+2));
