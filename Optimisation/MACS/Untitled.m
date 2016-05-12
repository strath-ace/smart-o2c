close all
clear all
clc


n_lambda = 50;
lambda = rand(100*n_lambda,3);

norms = (lambda(:,1).*lambda(:,1)+lambda(:,2).*lambda(:,2)+lambda(:,3).*lambda(:,3)).^0.5;

lambda = lambda./repmat(norms,1,3);

mins = min(lambda);
maxs = max(lambda);

[mm,ddt,eett,ene2t] = arch_shrk6([],[],lambda(1:n_lambda,:),0,[],mins,maxs,0,3,n_lambda); % first pass, add extremas
[mm,~,~,~] = arch_shrk6(mm,ddt,lambda(n_lambda+1:end,:),eett,ene2t,mins,maxs,0,3,n_lambda);

figure(1)

quiver3(zeros(size(lambda,1),1),zeros(size(lambda,1),1),zeros(size(lambda,1),1),lambda(:,1),lambda(:,2),lambda(:,3));
axis equal

figure(2)

quiver3(zeros(n_lambda,1),zeros(n_lambda,1),zeros(n_lambda,1),mm(:,1),mm(:,2),mm(:,3))
axis equal