function [memories,x,f,archsize]=upd_arch(memories,archsize,x,f,popsize,id,xsamp,fsamp,nsamp,lx,mfit)

for j=1:nsamp
    sempdom=archsize*ones(archsize,1);
    for k=1:archsize
        sempdom(k)=sum(fsamp(j,:)>= memories(k,lx+1:lx+mfit));
        if sempdom(k)==mfit
            break
        elseif sempdom(k)==0
            memories(k,lx+mfit+1)=memories(k,lx+mfit+1)+1;
        end
    end
    if max(sempdom)<mfit
        archsize=archsize+1;
        memories(archsize,1:lx)=xsamp(j,1:lx);
        memories(archsize,lx+1:lx+mfit)=fsamp(j,:);
        memories(archsize,1+lx+mfit)=0;
    end
    lock=0;
    for k=1:popsize
        if k~=id
            sempdom=sum(fsamp(j,:)>= f(k,:));
            if sempdom==0&&lock==0
                lock=1;
                f(k,:)=fsamp(j,:);
                x(k,:)=xsamp(j,:);
                break
            end
        end
    end
end

return