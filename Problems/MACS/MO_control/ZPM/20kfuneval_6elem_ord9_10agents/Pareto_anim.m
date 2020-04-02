% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
close all
clear
clc

load('memory.mat')


nit = max(mem.history(:,1));

mino1 = min(mem.history(:,end-3));
maxo1 = max(mem.history(:,end-3));
mino2 = min(mem.history(:,end-2));
maxo2 = max(mem.history(:,end-2));
deltao1 = maxo1-mino1;
deltao2 = maxo2-mino2;

figure()
box on
set(gca,'LooseInset',get(gca,'TightInset'));
set(gcf,'color','w');
%set(gca,'nextplot','replacechildren','visible','off');
plot(mem.history(mem.history(:,1)==i,end-3),mem.history(mem.history(:,1)==i,end-2),'b.','MarkerSize',10)
axis([mino1-deltao1*0.2 maxo1+deltao1*0.2 mino2-deltao2*0.2 maxo2+deltao2*0.2])
xlabel('-\theta_f [rad]')
ylabel('q [btu/ft^2/s]')

%[im,map] = rgb2ind(f.cdata,256,'nodither');
%im(1,1,1,40) = 0;
filename = 'Pareto_anim.gif';

for i = 1:nit
    
    plot(mem.history(mem.history(:,1)==i,end-3),mem.history(mem.history(:,1)==i,end-2),'b.','MarkerSize',10)
    str = ['iter ',num2str(i)];
    text(-0.3,70,str);
    axis([mino1-deltao1*0.2 maxo1+deltao1*0.2 mino2-deltao2*0.2 maxo2+deltao2*0.2])
    xlabel('10^{-6}\int u^Tu dt [lbf^2*ft^4/s^3]')
    ylabel('|h|_{max} [lbf*ft^2/s]')
    
    drawnow
    ax = gca;
    ax.Units = 'pixels';
    pos = ax.Position;
    ti = ax.TightInset;
        
    rect = [-ti(1), -ti(2), pos(3)+ti(1)+ti(3), pos(4)+ti(2)+ti(4)];
    f = getframe(ax,rect);
    
    im = frame2im(f);
    [imind,cm] = rgb2ind(im,256);
    
    if i == 1
        
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
        
    else
        
        imwrite(imind,cm,filename,'gif','WriteMode','append');
        
    end
    
end

plot(mem.history(mem.history(:,1)==i,end-3),mem.history(mem.history(:,1)==i,end-2),'b.','MarkerSize',10)
axis([mino1-deltao1*0.2 maxo1+deltao1*0.2 mino2-deltao2*0.2 maxo2+deltao2*0.2])

[~,b] = sort(mem(1).memory(:,end-2));
zz = mem(1).memory(b,:);

for i = 1:10
text(zz(i,end-3)+0.005,zz(i,end-2)+1,num2str(i));
end
str = ['iter ',num2str(40)];
text(-0.3,70,str);

xlabel('10^{-6}\int u^Tu dt [lbf^2*ft^4/s^3]')
ylabel('|h|_{max} [lbf*ft^2/s]')

drawnow
ax = gca;
ax.Units = 'pixels';
pos = ax.Position;
ti = ax.TightInset;

rect = [-ti(1), -ti(2), pos(3)+ti(1)+ti(3), pos(4)+ti(2)+ti(4)];
f = getframe(ax,rect);
im = frame2im(f);
[imind,cm] = rgb2ind(im,256);

   
imwrite(imind,cm,filename,'gif','WriteMode','append');

