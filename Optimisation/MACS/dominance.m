<<<<<<< HEAD:Optimisation/MACS/dominance.m
function dom=dominance(f,h)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Massimiliano Vasile ---------------------
%
% Computes dominance

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

[nf,mf]=size(f);

if (mf==1)                  %if problem is scalar
    [f,index]=sort(f);      %sort f
    dom(index)=0:1:nf-1;    %and dom is simply the order
else
    dom=zeros(nf,1);        %initialize a dominance vector
    if h<2                  

        for i=1:nf-1        %basically, for all couples (i,j)
            for j=i+1:nf
                %                if (i~=j)
                if all(f(i,:)==f(j,:))
                    if dom(i)<=dom(j)
                        dom(j)=nf*mf;
                    else
                        dom(i)=nf*mf;
                    end
                elseif (h==0)
                    %                        nz=find(f(i,:)~=f(j,:));
                    %                        if isempty(nz)==0
                    %                            if (all(f(i,nz)>f(j,nz)))
                    %[i,j f(i,:),f(j,:)]
                    %f(i,:)==f(j,:)
                    dom(i)=dom(i)+min(f(i,:)>=f(j,:));
                    dom(j)=dom(j)+min(f(j,:)>=f(i,:));
                    %dom(i)
                    %                                [dom(i),dom(j)]
                    %                                pause
                    %                            end
                    %                        end
                elseif (h==1)
                    %                        if (any(f(i,:)>f(j,:)))
                    dom(i)=dom(i)+sum(f(i,:)>f(j,:));
                    dom(j)=dom(j)+sum(f(j,:)>f(i,:));
                    
                    %                        end
                end

            end
        end

        %    end

    else
        dom=sum(f(1,:)==f(2,:))*max(f(1,:)>f(2,:))+sum(f(1,:)>f(2,:));
    end
end

return
=======
function dom=dominance(f,h)
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%-----------Copyright (C) 2016 University of Strathclyde-------------
%
%
%
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
% criterias. Moreover, the implemented way of computing the dominance can
% worsen the dominance index of two solutions with the same objective

[nf,mf]=size(f);

if (mf==1)                  %if problem is scalar
    [f,index]=sort(f);      %sort f
    dom(index)=0:1:nf-1;    %and dom is simply the order
else
    dom=zeros(nf,1);        %initialize a dominance vector
    if h<2                  

        for i=1:nf-1        %basically, for all couples (i,j)
            for j=i+1:nf
                %                if (i~=j)
                if all(f(i,:)==f(j,:))
                    if dom(i)<=dom(j)
                        dom(j)=nf*mf;
                    else
                        dom(i)=nf*mf;
                    end
                elseif (h==0)
                    %                        nz=find(f(i,:)~=f(j,:));
                    %                        if isempty(nz)==0
                    %                            if (all(f(i,nz)>f(j,nz)))
                    %[i,j f(i,:),f(j,:)]
                    %f(i,:)==f(j,:)
                    dom(i)=dom(i)+min(f(i,:)>=f(j,:));
                    dom(j)=dom(j)+min(f(j,:)>=f(i,:));
                    %dom(i)
                    %                                [dom(i),dom(j)]
                    %                                pause
                    %                            end
                    %                        end
                elseif (h==1)
                    %                        if (any(f(i,:)>f(j,:)))
                    dom(i)=dom(i)+sum(f(i,:)>f(j,:));
                    dom(j)=dom(j)+sum(f(j,:)>f(i,:));
                    
                    %                        end
                end

            end
        end

        %    end

    else
        dom=sum(f(1,:)==f(2,:))*max(f(1,:)>f(2,:))+sum(f(1,:)>f(2,:));
    end
end

return
>>>>>>> 5b7361d93c9119cf1d2e9e6c885bed93f924d71b:Optimisation/MACS/dominance.m
