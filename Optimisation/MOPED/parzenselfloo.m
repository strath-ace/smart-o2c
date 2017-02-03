% function y=parzenself(wei,x,f,ny,pdftype,rewei,tsp)
% wei: nx or empty row vector
% x: nx*mx matrix
% f: nx*mf matrix
% ny: scalar
% y: ny*mx matrix
% pdftype: string
% rewei: scalar

function y=parzenselfloo(wei,x,f,ny,pdftype,rewei,tsp)

if (nargin==0)
  disp('function y=parzenself(wei,x,f,ny,pdftype,rewei,tsp)');
  disp('wei: nx or empty row vector');
  disp('x: nx*mx matrix');
  disp('f: nx*mf matrix');
  disp('ny: scalar');
  disp('y: ny*mx matrix');
  disp('pdftype: string');
  disp('rewei: scalar');
  disp('tsp: scalar');
error(' ');
end;

[nx,mx]=size(x);

if size(f,1)==nx
else
  error('x and f must have the same number of rows');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wei_prova=wei;
%tsp=fix((rand*super+infer)*nx)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mf=size(f,2);
ind=zeros(ny,1);
y=zeros(ny,mx);

switch pdftype
  case 'norm'
    if rewei==0 
      if isempty(wei)
        kervar=cov(x,1)*nx^(-0.5/mx);
        ind=ceil(rand(ny,1)*nx);
      else
        swei=sum(wei);
        fmom=wei*x/swei;
        kervar=((repmat(wei,mx,1).*x')*x/swei-fmom'*fmom)*swei^(-0.5/mx);
        invvar=inv(kervar);
        zz=zeros(1,size(invvar,1));
        [scale,cost] = fminunc('normcost',zz,[],x,wei',invvar);
        if length(scale)==1
            scale=repmat(scale,length(invvar),1);
        end;
        if scale(1,1)==inf |scale(1,1)==nan
            'inf or nan'
        end
        scale=diag(sqrt(exp(scale)));
        invvar=scale*invvar*scale;
        kervar=inv(invvar);
        wei=cumsum(wei);
        wei=wei/wei(end);
        for i=1:ny
          ind(i)=ksearch(rand,wei);
        end;
      end;
    else
      if isempty(wei)
        sweif=nx;
        kervar=cov(f,1)*nx^(-0.5/mf);
        invwei=normparz(f,f,inv(kervar),repmat(1,1,nx));
      else
        sweif=sum(wei);
        fmom=wei*f/sweif;
        kervar=((repmat(wei,mf,1).*f')*f/sweif-fmom'*fmom)*sweif^(-0.5/mf);
        invwei=normparz(f,f,inv(kervar),wei);
      end;
      wei=(1./invwei').^rewei;
      swei=sum(wei);
      fmom=wei*x/swei;
      kervar=((repmat(wei,mx,1).*x')*x/swei-fmom'*fmom)*sweif^(-0.5/mx);
      wei=cumsum(wei);
      wei=wei/wei(end);
      for i=1:ny
        ind(i)=ksearch(rand,wei);
      end;
    end;
    [T,p]=chol(kervar);
    if p>0
        if max(kervar)==inf |kervar(1,1)==nan
            'inf'    
        end
      [T,p]=mvnfactor(kervar);
    end;
    if size(x,2)~=size(T,2)
        kervar=kervar-min(min(kervar))+1;
        [T,p]=chol(kervar);
        if p>0
            [T,p]=mvnfactor(kervar);
        end;
        'T<>x'
%        T
%        pause
%        x
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    ind=[];
%    for i_ind=1:ny
%        atsp=[];
%        cit=0;
%        while length(atsp)<tsp
%            pind=fix(rand.*nx)+1;
%            if atsp==[]
%                atsp=[atsp;pind];
%    else if sum(atsp==pind)<1
%                    atsp=[atsp;pind];
%                end
%            end
%            cit=cit+1;
%        end
%        wl=[];
%        for iw=1:tsp
%            wl=[wl;wei_prova(atsp(iw))];
%            if wei_prova(atsp(iw))==max(wl)
%                aind=atsp(iw);
%            end
%        end
%        ind=[ind;aind];
%    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%diag(kervar)'    
%ind
    y=randn(ny,size(T,1))*T+x(ind,:);
  case 't1'
    if rewei==0 
      if isempty(wei)
        kervar=var(x,1)*nx^(-0.5/mx);
        ind=ceil(rand(ny,1)*nx);
      else
        kervar=var(x,wei)*sum(wei)^(-0.5/mx);
        invvar=1./kervar;
        zz=zeros(1,size(invvar,1));
        [scale,cost] = fminunc('t1cost',0,[],x,invvar',wei');
        scale=exp(scale);
        if length(scale)==1
            scale=repmat(scale,length(invvar),1);
        end;
        invvar=scale'.*invvar;
        kervar=1./invvar;
        wei=cumsum(wei);
        wei=wei/wei(end);
        for i=1:ny
          ind(i)=ksearch(rand,wei);
        end;
      end;
    else
      if isempty(wei)
        sweif=nx;
        kervar=var(f,1)*nx^(-0.5/mf);
        invwei=t1parz(f,f,1./kervar,repmat(1,1,nx));
      else
        sweif=sum(wei);
        kervar=var(f,wei)*sweif^(-0.5/mf);
        invwei=t1parz(f,f,1./kervar,wei);
      end;
      wei=(1./invwei').^rewei;
      swei=sum(wei);
      kervar=var(x,wei)*sweif^(-0.5/mx);
      wei=cumsum(wei);
      wei=wei/wei(end);
      for i=1:ny
        ind(i)=ksearch(rand,wei);
      end;
    end;
    %if min(abs(kervar))<=0.0001
    %    kervar=kervar+.0001
    %end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    ind=[];
%    for i_ind=1:ny
%        atsp=[];
%        cit=0;
%        while length(atsp)<tsp
%            pind=fix(rand.*nx)+1;
%            if atsp==[]
%                atsp=[atsp;pind];
%    else if sum(atsp==pind)<1
%                    atsp=[atsp;pind];
%                end
%            end
%            cit=cit+1;
%        end
%        wl=[];
%        for iw=1:tsp
%            wl=[wl;wei_prova(atsp(iw))];
%            if wei_prova(atsp(iw))==max(wl)
%                aind=atsp(iw);
%            end
%        end
%        ind=[ind;aind];
%    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%kervar
%ind
    y=trnd(1,ny,mx)*diag(sqrt(kervar))+x(ind,:);
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function k=ksearch(refval,cumwei);

il=0;
iu=length(cumwei);

while (iu-il>1)
  im=floor((iu+il)/2);
  if (refval<cumwei(im))
    iu=im;
  else
    il=im;
  end;
end;

k=iu;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [T,p] = mvnfactor(sigma)
%MVNFACTOR  Do Cholesky-like decomposition, allowing zero eigenvalues
%   SIGMA must be symmetric.  In general T is not square or triangular.
%   P is the number of negative eigenvalues, and T is empty if P>0.

[U,D] = eig((sigma+sigma')/2);
D = diag(D);

tol = max(D) * length(D) * eps;
t = (abs(D) > tol);
D = D(t);
p = sum(D<0);

if (p==0)
   T = diag(sqrt(D)) * U(:,t)';
else
   T = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%