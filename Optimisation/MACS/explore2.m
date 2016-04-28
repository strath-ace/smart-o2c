function [xtrial,vtrial,ftrial,maxC,nfeval,discarded,rho]=explore2(compop,memories,x,v,f,cid,nfeval,lx,F,CR,rho,vlb,vub,coord_ratio,contr_ratio,rhoini,func,cp,cpat,ncon,arg)

% OLD INTERFACE IN HELP, BAD SIGN :-(
%
%  [x,v,f,nfeval]=explore(memories,box,idsamp_number,x,v,f,nfeval,lx,F,CR,archsize,popsize,vlb,vub,func,arg)
%
%
%  INPUT
%           memories       : archive
%           box            : boxes
%           idsamp_number  : selected box id
%           x              : positions
%           v              : velocities
%           f              : objective values
%           nfeval         : number of function evaluations
%           lx             : problem dimensions
%           F              : perturbation factor
%           CR             :
%           archsize       : archive size
%           popsize        : population size
%           vlb,vub        : boundaries of the search space
%           func           : function handle to objective function
%           arg            : arguments of func
%
%  OUTPUT
%           x       : updated x
%           v       : updated velocity
%           f       : updated function value
%           nfeval  : updated number of function evaluations

% CHANGE LOG
% Created by: Massimiliano Vasile 2010
% Revised, cleaned and optimized by: Lorenzo A. Ricciardi 2015


%% INITIALIZATION

vtrial=zeros(size(x));                                                      % vector of velocities of trial agents
ctrial=zeros(size(x,1),ncon);                                               % vector of contstraints of trial agents
maxC=zeros(size(x,1),1);                                                    % vector containing the max constraint violation of trial agents
xtrial=x;                                                                   % vector of position of trial agents, cloned from current position
ftrial=f;                                                                   % vector of objective function of trials agents, cloned from current one
discarded(length(compop)).f=[];                                             % structure containing objective function values of "discarded" agents
discarded(length(compop)).x=[];                                             % structure containing position of "discarded" agents
discarded(length(compop)).c=[];                                             % structure containing constraint violation of "discarded" agents
Delta=(vub-vlb)/2;                                                          % half span of search space
[npop,nf]=size(f);                                                          % number of agents within the population and dimensionality of criteria space (dimension of objective function f)
n_ids=ceil(lx*coord_ratio);                                                 % n_ids is equal to the number of dimensions which will actually be scanned (not automatically all of them, if coord_ratio <1 )


%% MAIN LOOP

for i=compop                                                                % for all elements in "compop" (need to check what does this mean)
    
    %% LOOP VARIABLES RESET
    
    impr=0;                                                                 % initialize impr to 0. It's a flag that will become 1 once an improvement has been made, and will prevent further search for this agent
    ids=randperm(lx);                                                       % create a vector from 1 to n_ids, and shuffle it's elements
    r=-1+2*rand(1,n_ids);                                                   % create a vector of random numbers between -1 and 1, with n_ids elements (so as the num of coordinates which will be employed in the local search)
    
    %% INERTIA (continue along the direction of the previous improvement)
    
    if norm(v(i,:))~=0                                                      % if velocity of that agent is not null (velocity is supplied as a parameter, so has been computed in previous steps). THIS IS ALWAYS FALSE DURING THE FIRST CALL OF EXPLORE
        
        alph=rand;                                                          % create another random number alpha between 0 and 1

        alph=alpha_clip(x(i,:),alph,v(i,:),vlb,vub);                        % new, faster way of clipping_alpha (still doesn't consider alternative behaviour, can be modified inside function)
        %[alph,v(i,:)]=alpha_clip2(x(i,:),alph,v(i,:),vlb,vub);             % modified behaviour function
        
        xtrial(i,:)=alph*v(i,:)+x(i,:);
        
        % evaluates move

        if norm(xtrial(i,:)-x(i,:))>0
            
            [discarded(i),impr] = trymove (xtrial(i,:),f(i,:),discarded(i),cp,cid(i),func,arg);
            nfeval=nfeval+1;
            
        end
    end
     
    %% DIFFERENTIAL (take the direction suggested by other agents)
    
    if impr==0
        
        ind = randperm(4);                                                  % index pointer array
        rot = (0:1:npop-1);                                                 % rotating index array (size NP)
        
        a1  = randperm(npop);                                               % shuffle locations of vectors
        rt = rem(rot+ind(1),npop);                                          % rotate indices by ind(1) positions
        a2  = a1(rt+1);                                                     % rotate vector locations
        rt = rem(rot+ind(2),npop);
        a3  = a2(rt+1);
        
        e = rand(1,lx) < CR;                                                % all random numbers < 0.5 are 1, 0 otherwise
        
        strategy='best';                                                    % this line should be a parameter supplied externally!!!
        
        switch strategy
        
            case 'best'
            
                af=randperm(nf);
                
                for j=1:nf
                
                    %                    [dummy,idf(j)]=min(f(:,j));
                    [~,idf(j)]=min(memories(:,lx+j));
                    
                end
                
                %xbest=x(idf(af(1)),:)
                xbest=memories(idf(af(1)),1:lx);
                dx =((xbest-x(i,:)) + F*(x(a1(i),:) - x(a2(i),:)));         % differential variation
                
            case 'rand'
                
                dx =(x(a3(i),:)-x(i,:) + F*(x(a1(i),:) - x(a2(i),:)));      % differential variation
                
        end
        
        xtrial(i,:)=e.*dx+x(i,:);
        
        % Double check feasibility error

        mask = (xtrial(i,:)<vlb)+(xtrial(i,:)>vub);                         % mask vector containing components out of bounds
        xnew = mask.*rand(1,lx).*(vub-vlb)+vlb;                             % vector of new components, where xtrial is out of bounds
        xtrial(i,:) = xtrial(i,:).*~mask+xnew.*mask;                        % new vector, completely in bounds

        % idea: why not repeat alpha clipping: we have a velocity vector (dx), and can consider alpha to be 1...
        
        % evaluates move

        if norm(xtrial(i,:)-x(i,:))>0

            [discarded(i),impr] = trymove (xtrial(i,:),f(i,:),discarded(i),cp,cid(i),func,arg);            
            nfeval=nfeval+1;
            
        else

            ftrial(i,:)=f(i,:);

        end

    end
    
    %% PATTERN SEARCH
    
    if impr==0
        
        xtrial(i,:)=cpat*xtrial(i,:)+(1-cpat)*x(i,:);                       % if cpat=1, start from xtrial coming from DE step, otherwise start from initial x (the one supplied as input parameter)
        %         ftrial(i,:)=f(i,:);
        %xsample=[];
        %fsample=[];
        %x0trial(i,:)=xtrial(i,:);
        %f0trial(i,:)=cpat*ftrial(i,:)+(1-cpat)*f(i,:);
        
        for j=1:n_ids                                                       % for every dimension (remember that they are shuffled devery time)
            
            xtrial(i,ids(j))=xtrial(i,ids(j))+r(j)*rho(i,1)*Delta(ids(j));  % take a step in that direction, with random amplitude between 0 and current max allowed by rho (either 1/2 of total span, 1/4 or 1/8)
            
            if xtrial(i,ids(j))<vlb(ids(j))                                 % check that the step is in the feasible area
            
                xtrial(i,ids(j))=vlb(ids(j));
                
            elseif xtrial(i,ids(j))>vub(ids(j))
                
                xtrial(i,ids(j))=vub(ids(j));
                
            end
            
            if norm(xtrial(i,:)-x(i,:))>0                                   % if agent has moved from previous position
                
                [discarded(i),impr] = trymove (xtrial(i,:),f(i,:),discarded(i),cp,cid(i),func,arg);                               
                nfeval=nfeval+1;
                
            end
                       
            if (all(ftrial(i,:)>=f(i,:))&&maxC(i)<=0)||(maxC(i)>0&&maxC(i)>cid(i))  % if previously sampled solution is worse than the previous (either for objectives or for constrains)
                
                %xsample=[xsample;xtrial(i,:) r(j)*rho(i,1)*Delta(ids(j))];  % xsample and fsample are used in the DDS phase, if needed
                %fsample=[fsample;ftrial(i,:)];                             
                
                if r(j)>0                                                   % if previous rand number was positive
                    
                    rr=rand;                                                % generate a new positive one
                    
                else                                                        % otherwise
                    
                    rr=-rand;                                               % generate a new negative one
                    
                end
                
                xtrial(i,ids(j))=x(i,ids(j))-rr*rho(i,1)*Delta(ids(j));     % the NEW trial position is on the OPPOSITE side of the previously sampled one (notice +r before and -rr now)
                
                if xtrial(i,ids(j))<vlb(ids(j))                             % always check bounds
                    
                    xtrial(i,ids(j))=vlb(ids(j));
                    
                elseif xtrial(i,ids(j))>vub(ids(j))
                    
                    xtrial(i,ids(j))=vub(ids(j));
                    
                end
                
                if norm(xtrial(i,:)-x(i,:))>0                               % if agent has moved from initial position
                    
                    [discarded(i),impr] = trymove (xtrial(i,:),f(i,:),discarded(i),cp,cid(i),func,arg);                   
                    nfeval=nfeval+1;
                    
                    % THIS SECTION IS VERY CONFUSED...
                    
                    if (all(ftrial(i,:)>=f(i,:))&&maxC(i)<=0)||(maxC(i)>0&&maxC(i)>cid(i))  % if still no improvement has been made
                        
                        %xsample=[xsample;xtrial(i,:) -rr*rho(i,1)*Delta(ids(j))];   % add previous sampled point to xsample
                        %fsample=[fsample;ftrial(i,:)];
                        
                        xtrial(i,ids(j))=x(i,ids(j));                       % reset xtrial to initial position
                        ftrial(i,:)=f(i,:);
                        
                    elseif (any(ftrial(i,:)<f(i,:))&&maxC(i)<=0)||(cid(i)>0&&maxC(i)<cid(i))    % if by some very mysterious reason (at this point) some improvement has been detected
                        
                        impr=1;
                        vtrial(i,:)=xtrial(i,:)-x(i,:);                     % update vtrial
                        
                    end
                    
                else                                                        % if no movement has been made
                    
                    xtrial(i,ids(j))=x(i,ids(j));                           % restore initial xtrial
                    ftrial(i,:)=f(i,:);
                    %xsample=[xsample;xtrial(i,:)  -rr*rho(i,1)*Delta(ids(j))];
                    %fsample=[fsample;ftrial(i,:)];
                    
                end
                
            elseif (any(ftrial(i,:)<f(i,:))&&maxC(i)<=0)||(cid(i)>0&&maxC(i)<cid(i))    % if an improvement was detected
                
                impr=1;
                vtrial(i,:)=xtrial(i,:)-x(i,:);                             % update vtrial
                
            end
            
            if impr==1                                                      % meh... totally illogic layout of code...
                
                break
                
            end
            
        end
        
    end
    
    %% DDS step
    
%     if impr==100&&any(i==id_pop_act_subpr)
%         
%         [vdds,J,nu]=dds(fsample,xsample,discarded(i).c,f0trial(i,:),x0trial(i,:),lambda(act_subpr(id_pop_act_subpr==i),:),1e-3,2);
%         tstep=rho(i,1);
%         delta_angle=0;
%         tol=1e-1;
%         step_trial=0;
%         while delta_angle<(1-tol)&&step_trial<10&&impr==0
%             
%             step_trial=step_trial+1;
%             dtrial=vdds;
%             for jj=1:lx
%                 if((dtrial(jj)+x0trial(i,jj))>vub(jj))&&(dtrial(jj)~=0)
%                     tstep=min([tstep,abs((vub(jj)-x0trial(i,jj))/dtrial(jj))]);
%                 elseif ((dtrial(jj)+x0trial(i,jj))<vlb(jj))&&(dtrial(jj)~=0)
%                     tstep=min([tstep,abs((x0trial(i,jj)-vlb(jj))/dtrial(jj))]);
%                 end
%             end
%             dtrial=vdds*tstep;
%             xtrial(i,:)=x0trial(i,:)+dtrial;
%             
%             if any(xtrial(i,:)<vlb)
%                 pause
%             end
%             if any(xtrial(i,:)>vub)
%                 pause
%             end
%             if cp==0
%                 ftrial(i,:)=func(xtrial(i,:),arg{:});
%                 maxC(i)=0;
%                 if ~all(ftrial(i,:)>=f(i,:))
%                     discarded(i).f=[discarded(i).f; ftrial(i,:)];
%                     discarded(i).x=[discarded(i).x; xtrial(i,:)];
%                     discarded(i).c=[discarded(i).c; maxC(i)];
%                     impr=1;
%                     
%                     vtrial(i,:)=xtrial(i,:)-x(i,:);
%                     
%                 end
%             else
%                 [ftrial(i,:),ctrial(i,:)]=func(xtrial(i,:),arg{:});
%                 maxC(i)=max(ctrial(i,:));
%                 if cid(i)<=0&&maxC(i)>0
%                     ftrial(i,:)=f(i,:)+maxC(i);
%                 end
%                 if (~all(ftrial(i,:)>=f(i,:))&&maxC(i)<=0)||(cid(i)>0&&maxC(i)<cid(i))
%                     discarded(i).f=[discarded(i).f; ftrial(i,:)];
%                     discarded(i).x=[discarded(i).x; xtrial(i,:)];
%                     discarded(i).c=[discarded(i).c; maxC(i)];
%                     impr=1;
%                     
%                     vtrial(i,:)=xtrial(i,:)-x(i,:);
%                     
%                 end
%             end
%             delta_angle=dot(ftrial(i,:)-f0trial(i,:),-lambda(act_subpr(id_pop_act_subpr==i),:));
%             delta_angle=delta_angle/(norm(ftrial(i,:)-f0trial(i,:))*norm(lambda(act_subpr(id_pop_act_subpr==i),:)));
%             tstep=tstep/2;
%             
%         end
%         
%     end
    
    %% MBH step_trial       WILL NEVER BE DONE WITH IMPRI==100!!!
    
%     if impr==100
%         %         kkk=0;
%         %         for kk=1:2:length(fsample(:,1))-1
%         %             kkk=kkk+1;
%         %             for l=1:mfit
%         %                 dxsamp(kkk,:,l)=quad1D([x(i,:);xsample(kk,:);xsample(kk+1,:)],[f(i,l);fsample(kk,l);fsample(kk+1,l)],0);
%         %             end
%         %         end
%         
%         %         [i pigr(act_subpr(id_pop_act_subpr==i))]
%         %         pause
%         if  any(i==id_pop_act_subpr)&&MBHflag>0&&pigr(act_subpr(id_pop_act_subpr==i))<0.8%<- TO BE REVISED WITH A PROPPER PARAMETER TO BE SET BY THE USER
%             niter=0;
%             fu=1e37;
%             
%             z=z-10;
%             fcurrent=g_MBHfun(x(i,:),func,lambda(act_subpr(id_pop_act_subpr==i),:),z,arg);
%             Z=f(i,:);
%             
%             while fu>=fcurrent&&niter<MBHflag
%                 niter=niter+1;
%                 foptionsNLP=optimset('Display','none','MaxFunEvals',100*lx,'LargeScale','off','Algorithm','interior-point','TolFun',1e-8);
%                 rvlb=x(i,:)-rho(i,1)*Delta;
%                 rvub=x(i,:)+rho(i,1)*Delta;
%                 rvlb=max([rvlb;vlb]);
%                 rvub=min([rvub;vub]);
%                 xtrial(i,:)=x(i,:)+(2*rand-1)*rho(i,1)*Delta;
%                 xtrial(i,:)=max([xtrial(i,:);rvlb]);
%                 xtrial(i,:)=min([xtrial(i,:);rvub]);
%                 
%                 ftrial(i,:)=func(xtrial(i,:),arg{:});
%                 
%                 Z=ftrial(i,:);
%                 xt=xtrial(i,:);
%                 %                 xt
%                 %                 c_MBHfun(xt,func,lambda(act_subpr(id_pop_act_subpr==i),:),Z,arg)
%                 %                  g_MBHfun(xt,func,lambda(act_subpr(id_pop_act_subpr==i),:),z,arg)
%                 % pause
%                 [u,fu,exitflag,output]=fmincon(@(xt)g_MBHfun(xt,func,lambda(act_subpr(id_pop_act_subpr==i),:),z,arg),xt,[],[],[],[],rvlb,rvub,@(xt)c_MBHfun(xt,func,lambda(act_subpr(id_pop_act_subpr==i),:),Z,arg),foptionsNLP);
%                 nfeval=nfeval+output.funcCount+1;
%                 if fu<fcurrent
%                     xtrial(i,:)=u;
%                     ftrial(i,:)=func(xtrial(i,:),arg{:});
%                     nfeval=nfeval+1;
%                     impr=1;
%                 end
%                 %                 lambda(act_subpr(id_pop_act_subpr==i),:)
%                 %                 fu-fcurrent
%                 %                 u
%                 %                 xtrial(i,:)
%                 %                 ftrial(i,:)
%                 %                 x(i,:)
%                 %                 f(i,:)
%                 %                 i
%                 %                g_fun(ftrial(i,:),lambda(act_subpr(id_pop_act_subpr==i),:),z)-                g_fun(f(i,:),lambda(act_subpr(id_pop_act_subpr==i),:),z)
%                 
%                 
%             end
%         end
%         
%         
%     end
    
    %% RADIUS CONTRACTION IN EVERYTHING FAILED, TO NARROW SEARCH AREA

    if impr==0                                                              % if nothing worked, contract rho
        
        rho(i,1)=rho(i,1)*contr_ratio;
        rho(i,2)=rho(i,2)+1;
        
        if rho(i,2)>3                                                       % unless it has contracted 3 times, in that case restart
            
            rho(i,1)=rhoini;
            rho(i,2)=0;
            
        end
        
    end
    
end


return
