function [J,grad,fobiet,gc]=calcola_Fc(xa,Y)

if isempty(Y.igen)
    igen=777;
else
    igen=Y.igen;
end
if isempty(Y.indexgrad)
    indexgrad=0;
else
    indexgrad=Y.indexgrad;
end
if isempty(Y.div)
    div=[];
else
    div=Y.div;
end


if indexgrad==0
    [m,n]=size(xa);
    [fobiet,gc]=calcola_F_Glen(xa,igen);
    J1(:,1)=fobiet(:,1);
    J1(:,2)=fobiet(:,2);
    for i=1:size(gc,2)
        J1(:,i+2)=gc(:,i);
    end
    
    grad=[];
    J=J1(1);
    CC=(gc(1,:)-div);
%    save constr CC
elseif indexgrad==1
    [m,n]=size(xa);
    delta=1e-3*ones(1,n);
    X=zeros(n,n);
    for i=1:n
        X(i,:)=xa(1,:);
        if X(i,i)<=1-delta(i)
            X(i,i)=X(i,i)+delta(i);
        else
            delta(i)=-1*delta(i);
            X(i,i)=X(i,i)+delta(i);
        end
    end
    Xc=[xa; X];
    [fobietX,gcX]=calcola_F_WindTurbine(Xc,igen);
    J1(1)=fobietX(1,1);
    J1(2)=fobietX(1,2);
    for i=1:size(gcX,2)
        J1(i+2)=gcX(1,i);
    end
    JX(:,1)=fobietX(2:end,1);
    JX(:,2)=fobietX(2:end,2);
    for i=1:size(gcX,2)
        JX(:,i+2)=gcX(2:end,i);
    end
    J=J1(1);
    grad=(JX(:,1)-J1(1))./delta';
    
elseif indexgrad==2
    [m,n]=size(xa);
    delta=1e-3*ones(1,n);
    X=zeros(n,n);
    for i=1:n
        X(i,:)=xa(1,:);
        if X(i,i)<=1-delta(i)
            X(i,i)=X(i,i)+delta(i);
        else
            delta(i)=-1*delta(i);
            X(i,i)=X(i,i)+delta(i);
        end
    end
    Xc=[xa; X];
    [fobietX,gcX]=calcola_F_WindTurbine(Xc,igen);
    J1(1)=fobietX(1,1);
    J1(2)=fobietX(1,2);
    for i=1:size(gcX,2)
        J1(i+2)=gcX(1,i);
    end
    JX(:,1)=fobietX(2:end,1);
    JX(:,2)=fobietX(2:end,2);
    for i=1:size(gcX,2)
        JX(:,i+2)=gcX(2:end,i);
    end
    J=J1(1);
    grad=(JX(:,1)-J1(1))./delta';
    CC=(gcX(1,:)-div);
    Ceq=[];
    CCX=(gcX(2:end,:)-repmat(div,n,1));
    gradCC=(CCX-repmat(CC,n,1))./repmat(delta',1,size(gcX,2));
    save constr CC gradCC
end





