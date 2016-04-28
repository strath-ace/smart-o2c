function [memories,dd,energy,ener2]=arch_rem(memories,dd,to_delete)

% Removes entries from archive, updating distances and energy
%
% memories=arch_rem(memories,dd,candidates,energy,lx,mfit,archsize,max_arch)
%
% INPUT
%       memories        :       archive
%       dd              :       matrix of pairwise energies
%       to_delete       :       agents to be removed from archive
%
% OUTPUT
%       memories        :       shrunk archive
%       dd              :       updated matrix of pairwise energies
%       energy          :       updated energy of the archive
%       ener2           :       new vector of energies of archive without
%                               i-th agent in archive
%
% Written by Lorenzo A. Ricciardi 2015

[memories,ids] = setdiff (memories,to_delete,'rows','stable');              % remove items from archive
dd = dd(ids,ids);                                                           % remove corresponding entries from dd matrix
energy = sum(sum(triu(dd,1)));                                              % update energy

ener2 = zeros(size(memories,1),1);

for i=1:size(ener2,1)
    
    ener2(i) = energy-sum(dd(i,1:end~=i));
    
end


end
