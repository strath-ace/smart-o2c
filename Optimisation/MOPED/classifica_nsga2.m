function [fobiet,x,fit,psc,funvin,cdis,GC]=classifica_nsga2(y,fob,fvin,cdisy,GCy,ncl,pam)

[Ni,Ng]=size(y);
Nf=size(fob,2);
cdis=cdisy;
GC=GCy;


ifvin0=find(fvin==0);
ifving0=find(fvin>0);
if ifvin0>=1
    con0=[y(ifvin0,:) fob(ifvin0,:) fvin(ifvin0)' ];
    f0 = non_domination_sort_mod(con0(:,1:end-1), Nf, Ng);
    [sf,isf]=sortrows([f0(:,end-1) -f0(:,end)]);
    pop=f0(isf,:);
else
    pop=[];
end

if ifving0>0
    cong0=[y(ifving0,:) fob(ifving0,:) fvin(ifving0)'];
    ncl=min(ncl,length(ifving0))
    ss=fix(size(cong0,1)/ncl);
    [sfvin,ifvin]=sort(fvin(ifving0)');
    for icl=1:ncl
        con(icl).F=cong0(ifvin((icl-1)*ss+1:icl*ss),:);
    end
    con(icl).F=[con(icl).F; cong0(ifvin(icl*ss+1:end),:)];
end

if ifving0>0
    for icl=1:ncl
        f0 = non_domination_sort_mod(con(icl).F(:,1:end-1), Nf, Ng);
        [sf,isf]=sortrows([f0(:,end-1) -f0(:,end)]);
        pop=[pop; f0(isf,:)];
    end
end

fobiet = pop(:,Ng+1:Ng+Nf);
x      = pop(:,1:Ng);
pam2 = 2 - pam;
fit = (pam2 - pam) / (Ni-1) * ((1:1:Ni)' - 1) + pam;
psc=pop(:,end-1);
try
for i=1:Ni
    for j=1:Ni
        if norm(y(i,:)-pop(j,1:Ng))<1e-16 && norm(fob(i,:)-pop(j,Ng+1:Ng+Nf))<1e-16
            funvin(j)=fvin(i);
            cdis(j,:)=cdisy(i,:);
            GC(j,:)=GCy(i,:);
        end
    end
end
catch
   disp('minkia') 
end




