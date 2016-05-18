function [J]=mgamsmask(y)
%global rplimit
%J = mgaCassini(y,rplimit);
 
        vlb =[-1000 30 100 30 400 1000]; 
        vub=[0 400 470 400 2000 6000]; 
y=(vub-vlb).*y+vlb;

[DVtot,rp] = mgaCASSINInoc(y);
w1=-sign(rp(1)-6351.8)+1;
w2=-sign(rp(2)-6351.8)+1;
w3=-sign(rp(3)-6778.1)+1;
w4=-sign(rp(4)-671492)+1;
 
J(1)=DVtot+0.01*(rp(1)-6351.8)^2*w1*0.5+0.01*(rp(2)-6351.8)^2*w2*0.5+0.01*(rp(3)-6778.1)^2*w3*0.5+0.001*(rp(4)-671492)^2*w4*0.5;
J(2)=sum(y(2:end));

return
