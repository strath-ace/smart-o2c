function [fobiet,funvin,funv,cdis]=constraints(pesi,div,fobiet,gc)

% Vincoli
funv=[];
cdis=[];
divn=div;
for j=1:size(fobiet,1)
    c=gc(j,:)-divn;
    cdis=[cdis;c];
    funv=[funv; max([zeros(1,length(c));c]).*pesi];
end


for j=1:size(fobiet,1)
    funv(j,:)=funv(j,:)./abs(divn);
    funvin(j)=sum(funv(j,:));
end