function memories=arch_shrk(memories,lx,mfit,max_arch)

% Archive shrinking to max allowed size
%
% memories=arch_shrk(memories,lx,mfit,archsize,max_arch)
%
% INPUT
%       memories        :       archive
%       lx              :       parameter space dimensions
%       mfit            :       objective function dimensions
%       max_arch        :       max allowed size
%
% OUTPUT
%       memories        :       shrunk archive
%
% Revised, cleaned and optimized by Lorenzo A. Ricciardi 2015

archsize = size(memories,1);                                                % number of agents currently in archive
[~,id_max]=min(memories(:,lx+1:lx+mfit),[],1);                              % compute ID of position of maximas for every objective value
id_max = unique(id_max);                                                    % beware of duplicates!

dd=Inf*ones(archsize);                                                      % matrix of square of distances (in criteria space) between agents

for ip=1:(archsize-1)                                                       % computation of matrix of square of distances (exploiting symmetry. Null distance is substituted with Inf)
    
    for jp=(ip+1):archsize
        
        dmem=((memories(ip,lx+1:lx+mfit)-memories(jp,lx+1:lx+mfit)));
        dd(ip,jp)=dmem*dmem';
        dd(jp,ip)=dd(ip,jp);
        
    end
    
end

selected=zeros(max_arch,1);                                                 % vector of chosen agents, with max_arch elements
selected(1:length(id_max))=id_max;                                          % the fist chosen agents are those which minimise the orthogonal subproblems
notselected=1:archsize;                                                     % vector of id of elements not selected, initially all current agents

notselected(selected(selected~=0)) = [];                                    % removing already chosen elements from list of notselected


for i=length(id_max)+1:max_arch
    
    [~,id2]=max(min(dd(selected(1:i-1),notselected)));                      % find the ID of the notselected agent that maximises the minimun distance from all selected agents
    selected(i)=notselected(id2);                                           % that agent is now selected
    notselected=notselected(notselected~=selected(i));                      % and removed from list of notselected agents
    
end

memories=memories(selected,:);

return