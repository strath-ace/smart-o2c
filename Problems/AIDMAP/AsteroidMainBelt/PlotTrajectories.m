function [r] = PlotTrajectories(Solutions,Costs,ListNodes)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%Addpath in case this file is ran seperately
addpath(genpath(fileparts(fileparts(fileparts(pwd)))));

dt = 1;
vvecsize = 3e7;
for i = 1:length(Solutions)
    r{i} = [];
    figure('Color',[1 1 1])
    axis equal; grid on; hold on
    for j = 2:length(Solutions{i})
        arrnode = ListNodes.(char(Solutions{i}(j)));
        depnode = ListNodes.(arrnode.parent);
        kep_trans = depnode.attributes.kep_trans;
        tvec = depnode.attributes.t_dep:dt:arrnode.attributes.t_dep;
        vvec = (arrnode.attributes.dV_dep.*vvecsize);
        [a,~] = size(r{i});
        for k = (1+a):(length(tvec)+a)
            r{i}(k,:) = StardustTool.CartesianElementsAt(kep_trans,tvec(k-a));
        end
        h = plot3(depnode.attributes.r_arr(1),depnode.attributes.r_arr(2),depnode.attributes.r_arr(3),'ko');
        set(h,'MarkerEdgeColor','none','MarkerFaceColor','k')
        
        quiver3(arrnode.attributes.r_dep(1),arrnode.attributes.r_dep(2),arrnode.attributes.r_dep(3),vvec(1),vvec(2),vvec(3),'g','linewidth',1.5,'MaxHeadSize',1);
        
%         if j == 2
%              plot3(arrnode.attributes.r_dep(1),arrnode.attributes.r_dep(2),arrnode.attributes.r_dep(3),'ko');
%         end
        
        temp = strsplit(depnode.node_ID,'____');
        temp = strsplit(char(temp(end)),'___');
        deptarget = temp(1);
        txt{i}(j) = text(depnode.attributes.r_arr(1),depnode.attributes.r_arr(2),depnode.attributes.r_arr(3),deptarget,'FontWeight','bold','VerticalAlignment','bottom');
        
        %dV Text
        dVtext = strcat(num2str(Costs{i}(j-1)),' km/s');
        txt2{i}(j) = text(arrnode.attributes.r_dep(1),arrnode.attributes.r_dep(2),arrnode.attributes.r_dep(3),dVtext,'FontWeight','bold','HorizontalAlignment','right');
        
        if j == length(Solutions{i})
            kep_trans = arrnode.attributes.kep_trans;
            tvec = arrnode.attributes.t_dep:dt:arrnode.attributes.t_arr;
            [a,~] = size(r{i});
            for k = (1+a):(length(tvec)+a)
                r{i}(k,:) = StardustTool.CartesianElementsAt(kep_trans,tvec(k-a));
            end
             h = plot3(arrnode.attributes.r_arr(1),arrnode.attributes.r_arr(2),arrnode.attributes.r_arr(3),'ko');
             %set('MarkerEdgeColor','k','MarkerFaceColor','none')
             set(h,'MarkerEdgeColor','k','MarkerFaceColor','k')
                          
             temp = strsplit(arrnode.node_ID,'____');
             temp = strsplit(char(temp(end)),'___');
             arrtarget = temp(1);
             
             txt{i+1}(j+1) = text(arrnode.attributes.r_arr(1),arrnode.attributes.r_arr(2),arrnode.attributes.r_arr(3),arrtarget,'FontWeight','bold','VerticalAlignment','bottom');
        end
            
    
    end
    plot3(r{i}(:,1),r{i}(:,2),r{i}(:,3),'b')
    xlabel('X [km]')
    ylabel('Y [km]')
    
    %Add total dV & ~ of asteroids to plot
    ylim=get(gca,'ylim');
    xlim=get(gca,'xlim');
    asteroidnum = length(Solutions{i})-1; %-1 for root
    dVtot = sum(Costs{i});
    txt3 = text(xlim(2),ylim(2),{char(strcat('Number of Asteroids:',{' '},num2str(asteroidnum))),char(strcat('Total \DeltaV:',{' '},num2str(dVtot),' km/s'))},'HorizontalAlignment','right','FontSize',16);
    

    
end

        
    
    


end

