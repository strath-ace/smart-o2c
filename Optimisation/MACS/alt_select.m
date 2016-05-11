close all
clear all
clc

n = 20;

x = rand(n,1);

ids = 1:n;
[xmin, idsxmin] = min(x);
[xmax, idsxmax] = max(x);

k = 7;

combins = combnk(ids,k);
combins = combins.*(repmat(((any(combins==idsxmin,2)).*(any(combins==idsxmax,2)))==1,1,k));

prm = combnk(1:k,2);

combins(combins(:,1)==0,:) = [];

nrg = zeros(size(combins,1),1);
length(nrg)


for p = 1:size(combins,1)
    nrg(p) = sum(1./abs(x(combins(p,prm(:,1)))-x(combins(p,prm(:,2)))));
end

[nrg,idold]=sort(nrg);
combins = combins(idold,:);

plot(x(combins(1,:)),ones(k,1).*nrg(1),'r*');
hold on
plot(x,ones(length(x))*nrg(1),'bo')

idsel = [idsxmin; idsxmax];
nrgsel= sum(1./abs(x(idsel(1))-x(idsel(2))));

idtry = ids;
idtry(idtry==idsxmin) = [];
idtry(idtry==idsxmax) = [];

for i=3:k
    nrgtry=zeros(length(idtry),1);
    for j=1:length(idtry)
        nrgtry(j)= sum(1./abs(x(idtry(j))-x(idsel)));
    end
    [nrgtry,id]=min(nrgtry);
    nrgsel = nrgsel+nrgtry;
    idsel=[idsel;idtry(id)];
    idtry(idtry==idtry(id)) = [];
end

xsel = x(idsel);
plot(xsel,ones(k,1).*nrgsel,'g^')
for i=1:k
    line([x(combins(1,i)); x(combins(1,i))], [nrg(1); nrgsel])
end

