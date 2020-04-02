% Downloaded from:
% http://web.mysites.ntu.edu.sg/epnsugan/PublicSite/Shared%20Documents/Forms/AllItems.aspx?RootFolder=%2fepnsugan%2fPublicSite%2fShared%20Documents%2fCEC%202011%2d%20RWP&FolderCTID=&View=%7bDAF31868%2d97D8%2d4779%2dAE49%2d9CEC4DC3F310%7d

function[]=display_plot(gbest,phi_desired,distance)

dim=length(gbest);
phi=linspace(0,360,5000);
yax(1)=array_factorcir(gbest,(pi/180)*phi(1),phi_desired,distance,dim);
maxi=yax(1);
for i=2:5000%This loop finds out the maximum gain 
    yax(i)=array_factorcir(gbest,(pi/180)*phi(i),phi_desired,distance,dim);
    if maxi<yax(i)
        maxi=yax(i);
    end;
end;
for i=1:5000%This loop normalizes the Y-axis and finds the normalized gain values in decibels 
    yax(i)=yax(i)/maxi;
    yax(i)=20*log10(yax(i));
end;
plot(phi,yax,'g')
xlabel('Azimuth angle(deg)');
ylabel('Gain(db)');