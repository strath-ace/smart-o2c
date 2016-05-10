function [memories,cnts]=arch_shrk4(memories,lx,mfit,max_arch)

% Archive shrinking to max allowed size, NSGA-II approach
%
% memories=arch_shrk4(memories,lx,mfit,archsize,max_arch)
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

dist = zeros(archsize,1);                                                   % vector of distance score

for i=1:size(f,2)                                                           % loop through all objectives
    
    [f_sort,idx] = unique(f(:,i));                                          % sorting of UNIQUE components
    
    dist(idx(1)) = Inf;                                                     % max score on extremes
    dist(idx(end)) = Inf;                                                   % max score on extremes
   
    f_delta = f_sort(end)-f_sort(1);                                        % normalization factor

    
    for i=2:length(f_sort)-1                                                % for all other agents
        
        dist(idx(i)) = dist(idx(i)) + ((f_sort(i+1)-f_sort(i-1))/f_delta);
        
    end
    
end

[dd,idx] = sort(dist,'descend');

idsel = idx(1:max_arch);

cnts = 0;
memories=memories(idsel,:);                                                 % finally, this is the shrunk archive

return