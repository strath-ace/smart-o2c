function [memories,nrg,cnts]=arch_shrk3(memories,lx,mfit,max_arch)

% Archive shrinking to max allowed size, with optimal disposition of agents
%
% memories=arch_shrk3(memories,lx,mfit,archsize,max_arch)
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
id_min = unique(id_min)';                                                   % beware of duplicates!

idsel = zeros(max_arch,1);                                                  % vector of IDs of selected points, initialized empty

xx = permute(f,[1,3,2]);                                                    % FAST computation of all relative square distances
xx = repmat(xx,1,size(f,1));
xxt = permute(xx,[2,1,3]);
dd = sum((xx-xxt).^2,3);

dd = dd.*(dd>0)+realmin.*(dd==0);                                           % Computation of pairwise energies, with clipping to 1/realmin=realmax for all elements with the same coordinates
dd = 1./dd;                                                                 

for i=1:length(id_min)                                                      % for sure, local extremas are included as the first elements, and never put in discussion
   
    idsel(i)=id_min(i); 
    
end

nrgsel= 0;                                                                  % computation of energy of the first two agents

if size(f,2)>1
    
    combins = nchoosek(id_min,2);                                               % compute all different combinations of minimas
    
    for i = 1:size(combins,1)
        
        nrgsel = nrgsel + dd(combins(i,1),combins(i,2));                        % add current couple energy
        
    end
    
end

idtry = ids;                                                                % id of tentative agents, initially all
idtry = setdiff(idtry,idsel);                                               % remove agents already selected (those which minimize orthogonal subproblems)
sums =0;
% first pass

for i=length(id_min)+1:max_arch                                             % for all remaining positions
    
    nrgtry = sum(dd(idsel(idsel~=0),idtry));                                % computation of all energies
    sums = sums+length(idsel(idsel~=0))*(length(idtry)-1);
    
    [nrgtry,id]=min(nrgtry);                                                % pick the one which minimizes total energy

    nrgsel = nrgsel+nrgtry;                                                 % update total energy
    idsel(i)=idtry(id);                                                     % add it to list of selected
    idtry(idtry==idtry(id)) = [];                                           % remove it from list of candidates
    
end

% refining passes
impr = 1;                                                                   % improved, initially 1, so refining will take place at least once
cnts = 0;

sums2 = 0;
while impr
    
    nrgsel_old = nrgsel;                                                    % reference to old energy, to understand if there's been an improvement
    
    for i=length(id_min)+1:max_arch                                         % for all remaining agents (i.e, all except minimas of single objectives)
    
        impr = 0;                                                           % set not improved
        nrgsel_tmp = nrgsel - sum(dd(idsel(i),idsel(idsel~=idsel(i)))) ;    % compute energy of system without this agent
        sums2 =  sums2 + length(idsel)-1;
        
        idtry = [idtry idsel(i)];                                           % re add this one to the list of possible agents
        idsel(i) = [];                                                      % remove it from the list of "good" ones
                
        nrgtry = sum(dd(idsel,idtry));                                      % computation of energies        
        
        sums2 =  sums2 + length(idsel(idsel~=0))*(length(idtry)-1);
        
        nrgtry = nrgsel_tmp+nrgtry;                                         % compute total energy of system in all these configurations
        
        [nrgtry,id]=min(nrgtry);                                            % pick the config of min energy
                
        nrgsel = min([nrgsel nrgtry]);                                      % the final energy is the min between the original one and the new min (i.e if there's an improvement pick the new one, otherwise the old one)
        idsel=[idsel;idtry(id)];                                            % add this agent to the list of good ones
        idtry(idtry==idtry(id)) = [];                                       % and remove it from list of candidates
        
    end
    
    impr = nrgsel<nrgsel_old;                                               % if there's been an improvement, tell it
    cnts = cnts+1;
    
end

memories=memories(idsel,:);                                                 % finally, this is the shrunk archive
nrg = nrgsel;                                                               % this is the energy

return