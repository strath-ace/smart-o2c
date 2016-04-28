function [v,J,lambda]=dds2(F,x,F0,x0,d,tol,type)

x = unique(x,'rows');
F = unique(F,'rows');
[nF,lF]=size(F);
%nF=nF/2;
% DF=zeros(lF,nF);
d=-d;
DF = [];
% J=DF;
lambda=d;
k=0;

for i=1:nF
    if type==2
        dx(i,:)=x(2+2*(i-1),1:end-1)-x(1+2*(i-1),1:end-1);
        %     k=find(dx(i,:));
        for j=1:lF
            DF(j,i)=(F(2+2*(i-1),j)-F(1+2*(i-1),j))/norm(dx(i,:));
        end
        dxn(i,:)=dx(i,:)/norm(dx(i,:));
    else
        if norm(x(i,1:end)-x0)>0
            k=k+1;
            dx(k,:)=x(i,1:end)-x0;
            %     k=find(dx(i,:));
            for j=1:lF
                DF(j,k)=(F(i,j)-F0(j))/norm(dx(k,:));
            end
            dxn(k,:)=dx(k,:)/norm(dx(k,:));
        end
    end
end

if ~isempty(DF) && all(all(~isnan(dxn)))
    J=DF*pinv(dxn');
    lambda=pinv(DF)*d';
    
    if cond(J)<tol
        [q,r]=qr(d);
        lambda=pinv(DF)*q(:,round(rand*(lF-2))+2);
    end
    
    v=lambda(1)*dxn(1,:);
    for i=2:length(lambda)
        v=v+lambda(i)*dxn(i,:);
    end
else
    
    v = zeros(size(x));
    J = 0;
    %lambda = 0;
end

% J=DF*pinv(dxn');
% lambda=pinv(DF)*d';
%
% if cond(J)<tol
%     [q,r]=qr(d);
%     lambda=pinv(DF)*q(:,round(rand*(lF-2))+2);
% end
%
% v=lambda(1)*dxn(1,:);
% for i=2:length(lambda)
%     v=v+lambda(i)*dxn(i,:);
% end

return