close all
clear all
clc

n = 500;

x = rand(n,1);
x = [x 1-log(x+1)];
%x = x./repmat(max(x),n,1);

ids = 1:n;
[xmin1, idsxmin1] = min(x(:,1));
[xmin2, idsxmin2] = min(x(:,2));

k = 50;

% tic
% combins = combnk(ids,k);
% combins = combins.*(repmat(((any(combins==idsxmin1,2)).*(any(combins==idsxmin2,2)))==1,1,k));
% 
% prm = combnk(1:k,2);
% 
% combins(combins(:,1)==0,:) = [];
% 
% nrg = zeros(size(combins,1),1);
% length(nrg)
% 
% 
% for p = 1:size(combins,1)
%     nrg(p) = sum(1./sum(abs(x(combins(p,prm(:,1)),:)-x(combins(p,prm(:,2)),:)),2));
% end
% toc
% 
% [nrg,idold]=sort(nrg);
% combins = combins(idold,:);
% 
% plot3(x(combins(1,:),1),x(combins(1,:),2),ones(k,1).*nrg(1),'r*');
% hold on

tic
idsel = zeros(k,1);
idsel(1) = idsxmin1; 
idsel(2) = idsxmin2;
nrgsel= sum(1./sum(abs(x(idsel(1),:)-x(idsel(2),:)),2));
plot(x(idsel(1:2),1),x(idsel(1:2),2),'bo','MarkerFaceColor',[1 0 0])
hold on
plot(x(:,1),x(:,2),'bo')
axis equal
drawnow
M = getframe;
idtry = ids;
idtry(idtry==idsxmin1) = [];
idtry(idtry==idsxmin2) = [];

for i=3:k
    nrgtry=zeros(length(idtry),1);
    for j=1:length(idtry)
        nrgtry(j) = sum(1./sum((repmat(x(idtry(j),:),size(idsel(idsel~=0),1),1)-x(idsel(idsel~=0),:)).^2,2));
        %nrgtry(j)= sum(1./sum(abs(repmat(x(idtry(j),:),size(idsel(idsel~=0),1),1)-x(idsel(idsel~=0),:)),2));
    end
    [nrgtry,id]=min(nrgtry);
    nrgsel = nrgsel+nrgtry;
    idsel(i)=idtry(id);
    idtry(idtry==idtry(id)) = [];
    %pause(0.25)
    plot(x(idsel(idsel~=0),1),x(idsel(idsel~=0),2),'bo','MarkerFaceColor',[1 0 0])
    hold on
    plot(x(:,1),x(:,2),'bo')
    axis equal
    drawnow
    M = [M; getframe];
end
toc

xsel = x(idsel,:);
%plot3(x(:,1),x(:,2),ones(size(x,1)).*2*nrgsel,'bo')

%plot3(xsel(:,1),xsel(:,2),ones(k,1).*nrgsel,'g^')
% plot(xsel(:,1),xsel(:,2),'*r')
% hold on
% plot(x(:,1),x(:,2),'bo')
% axis equal
% drawnow

% for i=1:k
%     line([x(combins(1,i),1); x(combins(1,i),1)], [x(combins(1,i),2); x(combins(1,i),2)], [nrg(1); nrgsel])
% end

% text(x(idsel(1),1),x(idsel(1),2),nrgsel,'0');
% text(x(idsel(2),1),x(idsel(2),2),nrgsel,'0');
% for i=3:length(idsel)
%     text(x(idsel(i),1),x(idsel(i),2),nrgsel,num2str(i-2))
% end

% next pass (only one)
refpass = 0;
impr = 1;

while impr
    %pause(0.25)
    refpass = refpass+1;
    nrgsel_old = nrgsel;
    for i=3:k
        impr = 0;
        nrgsel_tmp = nrgsel - sum(1./sum((repmat(x(idsel(i),:),k-1,1)-x(idsel(idsel~=idsel(i)),:)).^2,2));
        idtry = [idtry idsel(i)];   % re add this one to the list of possible agents
        idsel(i) = [];              % remove it from the list of "good" ones
        nrgtry=zeros(length(idtry),1);
        for j=1:length(idtry)
            nrgtry(j)= sum(1./sum((repmat(x(idtry(j),:),size(idsel(idsel~=0),1),1)-x(idsel(idsel~=0),:)).^2,2));
        end
        nrgtry = nrgsel_tmp+nrgtry;
        [nrgtry,id]=min(nrgtry);
        nrgsel = min([nrgsel nrgtry]);
        
        idsel=[idsel;idtry(id)];
        idtry(idtry==idtry(id)) = [];
    end
    impr = nrgsel<nrgsel_old;
    xsel= x(idsel,:);
    %plot3(xsel(:,1),xsel(:,2),ones(k,1).*nrgsel,'r*')
    hold off
    plot(xsel(:,1),xsel(:,2),'bo','MarkerFaceColor',[1 0 0])
    hold on
    plot(x(:,1),x(:,2),'bo')
    axis equal
    drawnow
    M = [M; getframe];
end

