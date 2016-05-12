close all
clear all
clc

load('global_3imp');
qq = qq(:,1:end-4);

qq = qq(:,2:end);

mins = min(qq(:,end-1:end),[],1);
maxs = max(qq(:,end-1:end),[],1);

[mem2,dd,energy,ener2,mins,maxs]=arch_shrk6([],[],qq(1:200,:),0,[],mins,maxs,5,2,200);

start = 201;
finished = 0;
while ~finished

    if start+10<=size(qq,1)
       
        endp = start+10;
        
    else
        
        endp=size(qq,1);
        finished = 1;
    end
    
    fprintf('%f %f\n',start,endp);
    
    [mem2,dd,energy,ener2,mins,maxs]=arch_shrk6(mem2,dd,qq(start:endp,:),energy,ener2,mins,maxs,5,2,200);

    start = start+10;
end