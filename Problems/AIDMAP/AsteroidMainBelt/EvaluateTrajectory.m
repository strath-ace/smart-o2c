function [r] = EvaluateTrajectory(Solutions,ListNodes,evalnum)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

dt = 1;
vvecsize = 5e6;   
p = 0;
for i = evalnum(1):evalnum(end)
    p = p + 1;
    r{p} = [];
    figure
    axis equal; grid on; hold on
    for j = 2:length(Solutions{i})
        arrnode = ListNodes.(char(Solutions{i}(j)));
        depnode = ListNodes.(arrnode.parent);
        kep_trans = depnode.attributes.kep_trans;
        tvec = depnode.attributes.t_dep:dt:arrnode.attributes.t_dep;
        vvec = (arrnode.attributes.dV_dep.*vvecsize);
        [a,b] = size(r{p});
        mink = 1+a;
        maxk = length(tvec)+a;
        for k = mink:maxk
            r{p}(k,:) = StardustTool.CartesianElementsAt(kep_trans,tvec(k-a));
        end
        plot3(depnode.attributes.r_arr(1),depnode.attributes.r_arr(2),depnode.attributes.r_arr(3),'rx');
        
        quiver3(arrnode.attributes.r_dep(1),arrnode.attributes.r_dep(2),arrnode.attributes.r_dep(3),vvec(1),vvec(2),vvec(3),'k','LineWidth',1.2);
        
        if j == 2
             plot3(arrnode.attributes.r_dep(1),arrnode.attributes.r_dep(2),arrnode.attributes.r_dep(3),'ko');
        end
        
        temp = strsplit(depnode.node_ID,'____');
        temp = strsplit(char(temp(end)),'___');
        deptarget = temp(1);
        txt{i}(j) = text(depnode.attributes.r_arr(1),depnode.attributes.r_arr(2),depnode.attributes.r_arr(3),deptarget);
        
        
        plot3(r{p}(mink:maxk,1),r{p}(mink:maxk,2),r{p}(mink:maxk,3))
        pause
        if j == length(Solutions{i})
            kep_trans = arrnode.attributes.kep_trans;
            tvec = arrnode.attributes.t_dep:dt:arrnode.attributes.t_arr;
            [a,b] = size(r{p});
            mink = 1+a;
            maxk = length(tvec)+a;
            for k = mink:maxk
                r{p}(k,:) = StardustTool.CartesianElementsAt(kep_trans,tvec(k-a));
            end
             plot3(arrnode.attributes.r_arr(1),arrnode.attributes.r_arr(2),arrnode.attributes.r_arr(3),'rx');
             
             temp = strsplit(arrnode.node_ID,'____');
             temp = strsplit(char(temp(end)),'___');
            arrtarget = temp(1);
             
             txt{i+1}(j+1) = text(arrnode.attributes.r_arr(1),arrnode.attributes.r_arr(2),arrnode.attributes.r_arr(3),arrtarget);
            
             plot3(r{p}(mink:maxk,1),r{p}(mink:maxk,2),r{p}(mink:maxk,3))
            pause
        end
            
        
    end
    
end

        
    
    


end

