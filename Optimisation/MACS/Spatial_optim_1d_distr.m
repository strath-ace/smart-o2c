close all
clear all
clc

n = 20;

x = rand(n,1);

ids = 1:n;
[xmin, idsxmin] = min(x);
[xmax, idsxmax] = max(x);

k = 4;

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

% combins = combins(round(linspace(1,length(nrg),100)),:);
% nrg = nrg(round(linspace(1,length(nrg),100)));

subplot(1,2,1)
for i = 1:length(nrg)
    plot(x(combins(i,:)),ones(k,1).*nrg(i),'bo');
    hold on
end

subplot(1,2,2)
plot(nrg)

% figure
% histogram(nrg)