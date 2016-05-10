close all
clear all
clc

n = 20;

x = rand(n,1);

ids = 1:n;
[xmin, idsxmin] = min(x);
[xmax, idsxmax] = max(x);

for k = 3:n-1
    
    combins = combnk(ids,k);
    combins = combins.*(repmat(((any(combins==idsxmin,2)).*(any(combins==idsxmax,2)))==1,1,k));
    
    prm = combnk(1:k,2);
    
    combins(combins(:,1)==0,:) = [];
    
    nrg = zeros(size(combins,1),1);
    
    
    for p = 1:size(combins,1)
        nrg(p) = sum(1./abs(x(combins(p,prm(:,1)))-x(combins(p,prm(:,2)))));
        
    end
    
    [~,id] = min(nrg);
    
    figure(k)
    xopt = x(combins(id,:));
    plot(x,zeros(n,1),'bo');
    hold on
    plot(xopt,zeros(k,1),'r*')
    
    figure(1000)
    nrg=sort(nrg);
    plot(nrg)
    hold on
    %axis([0 size(combins,1) 0 1000*nrg(1)])
end