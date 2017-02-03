function a=ckmin_2(a)
[u,d]=eig(a);
D=diag(d);

c=diag(a);
minval=min(abs(c));
    while sum(D<0)>=1
        c=diag(a);
        cmin=min(c);
        cmin2=(c==cmin);
        i=find(cmin2);
        a(i,i)=a(i,i)+minval*.01;
        [u,d]=eig(a);
        D=diag(d);
    end
