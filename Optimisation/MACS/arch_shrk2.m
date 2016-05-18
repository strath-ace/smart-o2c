function memories=arch_shrk2(memories,lx,mfit,max_arch)

% Archive shrinking to max allowed size, with optimal disposition of agents
%
% memories=arch_shrk2(memories,lx,mfit,archsize,max_arch)
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
% Written by Lorenzo A. Ricciardi 2015

archsize = size(memories,1);                                                % number of agents currently in archive
ids = 1:archsize;                                                           % ids of agents
f=memories(:,lx+1:lx+mfit);                                                 % objectives, easier to read

[~,id_min]=min(f,[],1);                                                     % compute ID of position of minimas for every objective value
id_min = unique(id_min)';                                                    % beware of duplicates!

idsel = zeros(max_arch,1);                                                  % vector of IDs of selected points, initialized empty

for i=1:length(id_min)                                                      % for sure, local extremas are included as the first elements, and never put in discussion
   
    idsel(i)=id_min(i); 
    
end

nrgsel= 1./sum((f(idsel(1),:)-f(idsel(2),:)).^2);                           % computation of energy of the first two agents

if length(id_min)>2                                                         % if more than 2 orthogonal subproblems, add energy of respective best agent
    
    for i=3:length(id_min)
    
        nrgsel = nrgsel + sum(1./sum((repmat(f(id_min(i),:),size(id_min(id_min~=id_min(i)),1),1)-f(id_min(id_min~=id_min(i)),:)).^2,2));
    
    end
    
end

idtry = ids;                                                                % id of tentative agents, initially all
idtry = setdiff(idtry,idsel);                                               % remove agents already selected (those which minimize orthogonal subproblems)

% first pass

for i=length(id_min)+1:max_arch                                             % for all remaining positions
        
    ftilde = permute(f(idtry,:),[3,2,1]);                                   % 3D matrix for maximum possible vectorization, beware, works but's really complex
    ftilde = repmat(ftilde,length(idsel(idsel~=0)),1,1);
    ftilde2 = f(idsel(idsel~=0),:);
    ftilde2 = repmat(ftilde2,1,1,size(ftilde,3));
    
    nrgtry = sum(1./sum((ftilde-ftilde2).^2,2));                            % computation of all energies
    
    [nrgtry,id]=min(nrgtry);                                                % pick the one which minimizes total energy

    nrgsel = nrgsel+nrgtry;                                                 % update total energy
    idsel(i)=idtry(id);                                                     % add it to list of selected
    idtry(idtry==idtry(id)) = [];                                           % remove it from list of candidates
    
end

% refining passes
impr = 1;                                                                   % improved, initially 1, so refining will take place at least once

while impr
    
    nrgsel_old = nrgsel;                                                    % reference to old energy, to understand if there's been an improvement
    
    for i=length(id_min)+1:max_arch                                         % for all remaining agents (i.e, all except minimas of single objectives)
    
        impr = 0;                                                           % set not improved
        nrgsel_tmp = nrgsel - sum(1./sum((repmat(f(idsel(i),:),max_arch-1,1)-f(idsel(idsel~=idsel(i)),:)).^2,2));  % compute energy of system without this agent
        idtry = [idtry idsel(i)];                                           % re add this one to the list of possible agents
        idsel(i) = [];                                                      % remove it from the list of "good" ones
        
        
        ftilde = permute(f(idtry,:),[3,2,1]);                               % 3D matrix for maximum possible vectorization, beware, works but's really complex
        ftilde = repmat(ftilde,length(idsel(idsel~=0)),1,1);
        ftilde2 = f(idsel(idsel~=0),:);
        ftilde2 = repmat(ftilde2,1,1,size(ftilde,3));
        
        nrgtry = sum(1./sum((ftilde-ftilde2).^2,2));                        % computation of energies
        
        
        nrgtry = nrgsel_tmp+nrgtry;                                         % compute total energy of system in all these configurations
        [nrgtry,id]=min(nrgtry);                                            % pick the config of min energy
        
        
        nrgsel = min([nrgsel nrgtry]);                                      % the final energy is the min between the original one and the new min (i.e if there's an improvement pick the new one, otherwise the old one)
        idsel=[idsel;idtry(id)];                                            % add this agent to the list of good ones
        idtry(idtry==idtry(id)) = [];                                       % and remove it from list of candidates
        
    end
    
    impr = nrgsel<nrgsel_old;                                               % if there's been an improvement, tell it
    
end

memories=memories(idsel,:);                                                 % finally, this is the shrunk archive

return