function [memories,dd,energy,ener2]=arch_shrk5(memories,dd,candidates,energy,ener2,lx,mfit,max_arch)

% Archive shrinking to max allowed size, with optimal disposition of agents
%
% memories=arch_shrk5(memories,dd,candidates,energy,lx,mfit,archsize,max_arch)
%
% INPUT
%       memories        :       archive
%       dd              :       matrix of pairwise energies
%       candidates      :       new agents to be added, if possible
%       energy          :       energy of the archive
%       lx              :       parameter space dimensions
%       mfit            :       objective function dimensions
%       max_arch        :       max allowed size
%
% OUTPUT
%       memories        :       shrunk archive
%       dd              :       updated matrix of pairwise energies
%       energy          :       updated energy of the archive
%
% Written by Lorenzo A. Ricciardi 2015

% Memories is assumed to be self dominated and not dominated by any
% candidate
% Candidates are assumed to be self dominated and not dominated by any
% memory
% This mean that the archive cannot decrease within this function
% This also means that candidates CANNOT be new minimas/maximas (otherwise
% there would be domination)

oldsize = size(memories,1);                                                 % number of agents currently in archive
candsize = size(candidates,1);                                              % number of candidates
qty_to_add = min(max_arch-oldsize,candsize);                                % take all candidates if they fit into the archive, or less if there's too many of them

if qty_to_add == candsize                                                   % all candidates will be archived, no selection performed, only update of energy and distance matrix
    
    f = candidates(:,lx+1:lx+mfit);                                         % objective values of candidates
    
    xx = permute(f,[1,3,2]);                                                % FAST computation of all relative square distances of candidates
    xx = repmat(xx,1,size(f,1));
    xxt = permute(xx,[2,1,3]);
    dd_new = sum((xx-xxt).^2,3);
    
    dd_new = dd_new.*(dd_new>0)+realmin.*(dd_new==0);                       % Computation of pairwise energies, bottom right diagonal block
    dd_new = 1./dd_new;
    
    new_energy_contr = sum(sum(triu(dd_new,1)));                            % Energy contribution of bottom right diagonal block matrix
    
    dold = dd;                                                              % copy of old dd
    dd = zeros(oldsize+qty_to_add,oldsize+qty_to_add);                      % new (larger) dd
    
    dd(1:oldsize,1:oldsize) = dold;                                         % inserting old dd into new one
    
    dd(oldsize+1:oldsize+qty_to_add,oldsize+1:oldsize+qty_to_add) = dd_new; % adding dd_new, bottom right block matrix
    energy = energy+new_energy_contr;                                       % update total energy
    
    if oldsize>0                                                            % if there were previous elements in archive
        
        fold = memories(:,lx+1:lx+mfit);                                    % objective values of agents in memory
        
        xx = permute(fold,[1,3,2]);                                         % FAST computation of pairwise energies, top right block matrix
        xx = repmat(xx,1,candsize);
        yy = permute(f,[3,1,2]);
        yy = repmat(yy,oldsize,1);
        drect = 1./sum((xx-yy).^2,3);
        
        dd(1:oldsize,oldsize+1:oldsize+qty_to_add) = drect;                 % adding bloc matrices to new dd
        dd(oldsize+1:oldsize+qty_to_add,1:oldsize) = drect';
        
        cross_ener_contr = sum(sum(drect));                                 % energy contribution of block top right matrix (cross contribution)
        
        energy = energy+cross_ener_contr;                                   % update energy
        
    end
    
    memories = [memories;candidates];                                       % store all candidates
    
    ener2 = zeros(size(memories,1),1);
    
    for i = 1:size(memories,1)
        
        ener2(i) = energy-sum(dd(i,1:end~=i));
        
    end
    
end

if (qty_to_add < candsize)                                                  % not all candidates will fit into the archive, selection will be performed
    
    remaining_candidates = 1:candsize;                                      % ids of candidates to be evaluated
    
    if qty_to_add >0
        
        memories = [memories; zeros(qty_to_add,size(memories,2))];              % enlarge memories
        
        dold = dd;                                                              % copy of old dd
        dd = zeros(oldsize+qty_to_add,oldsize+qty_to_add);                      % new (larger) dd
        
        dd(1:oldsize,1:oldsize) = dold;                                         % inserting old dd into new one
        
        % selection of best candidates while keeping current archive, i.e we
        % add candidates which minimize total energy growth
        
        for i = 1:qty_to_add                                                    % for all vacancies
            
            f = candidates(remaining_candidates,lx+1:lx+mfit);                  % objective values of remaining candidates
            fold = memories(1:oldsize,lx+1:lx+mfit);                           % objective values of agents in memory
            
            xx = permute(fold,[1,3,2]);                                         % FAST computation of pairwise energies, top right block matrix
            xx = repmat(xx,1,size(f,1));
            yy = permute(f,[3,1,2]);
            yy = repmat(yy,size(fold,1),1);
            drect = 1./sum((xx-yy).^2,3);
            
            [en_best_cand,id_best_cand] = min(sum(drect));                      % location of best candidate
            
            memories(oldsize+1,:) = candidates(remaining_candidates(id_best_cand),:); % add it to memories
            
            dd(1:oldsize,oldsize+1) = drect(:,id_best_cand);                    % adding relevant drect entry to dd
            dd(oldsize+1,1:oldsize) = drect(:,id_best_cand)';                   %
            
            dd(oldsize+1,oldsize+1) = 1/realmin;                                % bottom right element of dd (new one)
            
            energy = energy + en_best_cand;                                     % update energy
            
            remaining_candidates(id_best_cand) = [];                            % remove it from remaining candidates, to avoid chosing duplicates
            
            oldsize = oldsize+1;                                                % size of archive has just increased by one
            
        end
        
        for i = 1:size(memories,1)

            ener2(i) = energy-sum(dd(i,1:end~=i));
            
        end
        
    end
    

    if qty_to_add == 0
        
        cnt = 0;
        
        impr = 1;                                                               % set improved just to enter the loop
        
        while impr && cnt<100                                                   % if there's been an improvement over the last loop, try again
            
            impr = 0;                                                           % assume no improvement has been made (avoids infinite loop)
            id_best_cand = zeros(max_arch,1);                                   % vector of id of best candidate replacing ith agent of the archive
            en_best_cand = zeros(max_arch,1);                                   % vector of total energies associated to the replacement of ith agent with most promising candidate
            
            
            for i = 1:max_arch                                                  % for all elements in the archive
                
                tmp_nrg = ener2(i);                                             % energy of archive without this element
                
                f = candidates(remaining_candidates,lx+1:lx+mfit);              % objective values of remaining candidates
                fold = memories(1:end~=i,lx+1:lx+mfit);                         % objective values of agents in memory
                
                xx = permute(fold,[1,3,2]);                                     % FAST computation of pairwise energies, top right block matrix
                xx = repmat(xx,1,size(f,1));
                yy = permute(f,[3,1,2]);
                yy = repmat(yy,size(fold,1),1);
                drect = 1./sum((xx-yy).^2,3);
                
                [en_best_cand(i),id_best_cand(i)] = min(sum(drect)+tmp_nrg);   % location of best candidate
                
            end
            
            [en_best_global, id_best_global] = min(en_best_cand);               % location of best combination
            
            % thus
            % best result is substitute id_best_global with
            % id_best_cand(id_best_global)
            
            if en_best_global<energy
                
                impr = 1;
                
                tmp_mem = memories(id_best_global,:);                           % store agent, will be put into the candidates
                memories(id_best_global,:) = candidates(remaining_candidates(id_best_cand(id_best_global)),:);% swapping the two agents
                candidates(remaining_candidates(id_best_cand(id_best_global)),:) = tmp_mem; % put temp element into the candidates, for possible re-pick
                
                energy = en_best_global;                                        % updating energy
                
                f = memories(1:end~=id_best_global,lx+1:lx+mfit);               % FAST computation of energies
                thisf = memories(id_best_global,lx+1:lx+mfit);
                
                xx = permute(f,[1 3 2]);
                yy = permute(thisf,[1 3 2]);
                yy = repmat(yy,size(f,1),1);
                
                thisdd = 1./sum((xx-yy).^2,3);
                
                c=false(1,length(thisdd)+1);                                    % inserting 1/realmin in appropriate position
                c(id_best_global)=true;
                result=nan(size(c));
                result(~c)=thisdd;
                result(c)=1/realmin;
                
                dd(id_best_global,:) = result;                                  % update dd
                dd(:,id_best_global) = result';
                
                for i = 1:size(memories,1)
                    
                    ener2(i) = energy-sum(dd(i,1:end~=i));
                    
                end
                
            end
           
            cnt = cnt+1;
            
        end
        
    end
    
end

end