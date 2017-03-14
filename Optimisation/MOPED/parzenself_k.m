% function y=parzenself(wei,x,f,ny,pdftype,rewei)
% wei: nx or empty row vector
% x: nx*mx matrix
% f: nx*mf matrix
% ny: scalar
% y: ny*mx matrix
% pdftype: string
% rewei: scalar

function y=parzenself_k(wei,x,f,ny,pdftype,rewei)

if (nargin==0)
  disp('function y=parzenself(wei,x,f,ny,pdftype,rewei)');
  disp('wei: nx or empty row vector');
  disp('x: nx*mx matrix');
  disp('f: nx*mf matrix');
  disp('ny: scalar');
  disp('y: ny*mx matrix');
  disp('pdftype: string');
  disp('rewei: scalar');
  error(' ');
end;

[nx,mx]=size(x);

if size(f,1)==nx
else
  error('x and f must have the same number of rows');
end;

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
        %a=min(diag(kervar));
        %if a<1e-3
            %e_a=log10(1e-3)/log10(a);
            %b1=(kervar>=0);
            %b=b1.*kervar;
            %c=(kervar<0).*kervar;
            %kervar=b.^e_a-((-c).^e_a);
            %kervar=kervar./.001;
            %end
        % -_-_-_-_-_-_
        %valparz=normparz(x,x,inv(kervar),wei);
        %imwei=find(wei==max(wei));
        %imval=find(valparz==max(valparz));
        %while imwei~=imval
        %    wei=wei.^1.1;
        %    swei=sum(wei);
        %    fmom=wei*x/swei;
        %    kervar=((repmat(wei,mx,1).*x')*x/swei-fmom'*fmom)*swei^(-0.5/mx);
        %    valparz=normparz(x,x,inv(kervar),wei);
        %    imwei=find(wei==max(wei));
        %    imval=find(valparz==max(valparz));
        %end
        % -_-_-_-_-_-_
        diag(kervar)';
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
      [T,p]=mvnfactor(kervar);
    end;
    indice=0;
    while size(x,2)~=size(T,2) & indice<=20
        indice=indice+1
        kervar=kervar+max([1e-15;abs(min(min(kervar)))]);
        [T,p]=chol(kervar);
        if p>0
            [T,p]=mvnfactor(kervar);
        end;
        'T<>x'
%        T
%        pause
%        x
    end
    if indice>0
        kervar
    end
    y=randn(ny,size(T,1))*T+x(ind,:);
    for i_ind=1:size(y,1)
        for j_gen=1:size(y,2)    
            if (y(i_ind,j_gen)<0), y(i_ind,j_gen)= x(ind(i_ind),j_gen)/2; end,
            if (y(i_ind,j_gen)>1), y(i_ind,j_gen)=(x(ind(i_ind),j_gen)+1)/2; end,
        end
    end
  case 't1'
    if rewei==0 
      if isempty(wei)
        kervar=var(x,1)*nx^(-0.5/mx);
        ind=ceil(rand(ny,1)*nx);
      else
        kervar=var(x,wei)*sum(wei)^(-0.5/mx);
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
    y=trnd(1,ny,mx)*diag(sqrt(kervar))+x(ind,:);
      for i_ind=1:size(y,1)
        for j_gen=1:size(y,2)    
            if (y(i_ind,j_gen)<0), y(i_ind,j_gen)= x(ind(i_ind),j_gen)/2; end,
            if (y(i_ind,j_gen)>1), y(i_ind,j_gen)=(x(ind(i_ind),j_gen)+1)/2; end,
        end
    end
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