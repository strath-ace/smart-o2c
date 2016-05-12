function dom=dominance(f,h)
%
%   dominance
%
%  [dom]=dominance(f,h)
%
%  Evaluates the dominance of  element i
%  over any other element j in f
%
%   f(i,:)>f(j,:)
%
%   if h=0
%      simple dominance is computed, i.e.
%      i is dominated by j if all the not equal components of j are better than i
%   if h=1
%      strong dominance, i.e.
%      i is dominated by j if any component of j is better than i and the
%      dominance level is the sum of all the better components
%
%   (c)  Massimiliano Vasile  2003
% Modified: F.Zuiani 3/3/2010, changed behaviour for single objective
% functions. Now outputs cardinal ranking.
% Modified: F.Zuiani 12/3/2010, changed behaviour when two individuals i,j 
% are identical. Now one of them is not penalized while a very high
% dominance index is assigned to the other.
% Note: L.Ricciardi 15/01/2015, the penalization described above can preclude 
% diversity in parameter space and thus exploration capabilities, 
% since two different vectors in parameter space can lead to the same
% criterias. 

if h>=2
   error('Dominance should be either normal (0) or strong (1)') 
end

[nf,mf]=size(f);

if (mf==1)                                                                  % if problem is scalar
    [~,index]=sort(f);                                                      % sort f
    dom(index)=0:1:nf-1;                                                    % and dom is simply the order
else
    dom=zeros(nf,1);                                                        % initialize a dominance vector
    pen=zeros(nf,1);
end

if nf==1
    dom=0;
    return
end

if h==0
    for i=1:nf                                                              % for every agent
        ids = 1:nf;
        ids(i) = [];
        qq = repmat(f(i,:),nf-1,1);
        rr = f(ids,:);
        dom(ids) = dom(ids)+any(qq~=rr,2).*all(qq<=rr,2)+~pen(ids).*all(qq==rr,2).*nf.*mf;
        pen(ids) = pen(ids) + i.*all(qq==rr,2);
    end
else                                                     
    for i=1:nf                                                              % for every agent
        ids = 1:nf;
        ids(i) = [];
        qq = repmat(f(i,:),nf-1,1);
        rr = f(ids,:);
        dom(ids) = dom(ids)+any(qq~=rr,2).*all(qq<rr,2)+~pen(ids).*all(qq==rr,2).*nf.*mf;
        pen(ids) = pen(ids) + i.*all(qq==rr,2);
    end    
end

for i=1:nf
    if pen(i)~=0 && i==pen(pen(i))
        pen(i)=0;
        dom(i) = 0;
    end
end
    
return
