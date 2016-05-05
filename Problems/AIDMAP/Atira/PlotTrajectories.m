function [r] = PlotTrajectories(Solutions,ListNodes)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

dt = 1;
vvecsize = 5e6;
for i = 1:length(Solutions)
    r{i} = [];
    figure
    axis equal; grid on; hold on
    for j = 2:length(Solutions{i})
        arrnode = ListNodes.(char(Solutions{i}(j)));
        depnode = ListNodes.(arrnode.parent);
        kep_trans = depnode.attributes.kep_trans;
        tvec = depnode.attributes.t_dep:dt:arrnode.attributes.t_dep;
        vvec = (arrnode.attributes.dV_dep.*vvecsize);
        %vvec = (arrnode.attributes.dV_dep./arrnode.attributes.dV_tot)*vvecsize;
        [a,b] = size(r{i});
        for k = (1+a):(length(tvec)+a)
            r{i}(k,:) = StardustTool.CartesianElementsAt(kep_trans,tvec(k-a));
        end
        plot3(depnode.attributes.r_arr(1),depnode.attributes.r_arr(2),depnode.attributes.r_arr(3),'rx');
        if j ~= length(Solutions{i})
            quiver3(arrnode.attributes.r_dep(1),arrnode.attributes.r_dep(2),arrnode.attributes.r_dep(3),vvec(1),vvec(2),vvec(3),'g');
        end
        if j == 2
             plot3(arrnode.attributes.r_dep(1),arrnode.attributes.r_dep(2),arrnode.attributes.r_dep(3),'ko');
        end
        temp = strsplit(depnode.node_ID,'____');
        temp = strsplit(char(temp(end)),'___');
        deptarget = temp(1);
        txt{i}(j) = text(depnode.attributes.r_arr(1),depnode.attributes.r_arr(2),depnode.attributes.r_arr(3),deptarget);
    end
    plot3(r{i}(:,1),r{i}(:,2),r{i}(:,3),'b')
    
end

        
    
    


end

