function [xtrial,vtrial,ftrial,maxC,nfeval,discarded,rho]=exploreold(compop,memories,x,v,f,cid,...
         id_pop_act_subpr,act_subpr,lambda,z,nfeval,lx,mfit,F,CR,rho,vlb,vub,coord_ratio,contr_ratio,tolconv,rhoini,func,cp,pigr,MBHflag,cpat,ncon,arg)
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


vtrial=zeros(size(x));
ctrial=zeros(size(x,1),ncon);
maxC=zeros(size(x,1),1);
xtrial=x;
ftrial=f;
discarded.f=[];
discarded.x=[];
discarded.c=[];
Delta=(vub-vlb)/2;
[npop,nf]=size(f);

succmoves = zeros(size(x,1),1);
for i=compop
    nfbeg = nfeval;
    
    
    discarded(i).f=[];
    discarded(i).x=[];
    discarded(i).c=[];
    
    n_ids=ceil(lx*coord_ratio);
    impr=0;
    ids=randperm(lx);
    r=-1+2*rand(1,n_ids);
    %%
    %          INERTIA
    %          continue along the direction of the previous improvement
    if norm(v(i,:))~=0
        alph=rand;
        % Computes contraction of the step length
        for jj=1:lx
            if((v(i,jj)+x(i,jj))>vub(jj))&&(v(i,jj)~=0)
                alph=min([alph,abs((vub(jj)-x(i,jj))/v(i,jj))]);
            elseif ((v(i,jj)+x(i,jj))<vlb(jj))&&(v(i,jj)~=0)
                alph=min([alph,abs((x(i,jj)-vlb(jj))/v(i,jj))]);
            end
        end
        
        xtrial(i,:)=alph*v(i,:)+x(i,:);
        
        % Double check feasibility error
        for jj=1:lx
            if((xtrial(i,jj)-vlb(jj))<0)
                xtrial(i,jj)=vlb(jj);
            elseif((xtrial(i,jj)-vub(jj))>0)
                xtrial(i,jj)=vub(jj);
            end
        end
        % evaluates move
        if norm(xtrial(i,:)-x(i,:))>0
            
            if cp==0
                ftrial(i,:)=func(xtrial(i,:),arg{:});
                maxC(i)=0;
                if ~all(ftrial(i,:)>=f(i,:))
                    discarded(i).f=[discarded(i).f; ftrial(i,:)];
                    discarded(i).x=[discarded(i).x; xtrial(i,:)];
                    discarded(i).c=[discarded(i).c; maxC(i)];
                    impr=1;
                end
            else
                [ftrial(i,:),ctrial(i,:)]=func(xtrial(i,:),arg{:});
                maxC(i)=max(ctrial(i,:));
                if cid(i)<=0&&maxC(i)>0
                    ftrial(i,:)=f(i,:)+maxC(i);
                end
                if (~all(ftrial(i,:)>=f(i,:))&&maxC(i)<=0)||(cid(i)>0&&maxC(i)<cid(i))
                    discarded(i).f=[discarded(i).f; ftrial(i,:)];
                    discarded(i).x=[discarded(i).x; xtrial(i,:)];
                    discarded(i).c=[discarded(i).c; maxC(i)];
                    impr=1;
                end
            end
            
            nfeval=nfeval+1;
            
        end
    end
     
    %%
    %      DIFFERENTIAL
    if impr==0
        %    take the direction suggested by other agents
        ind = randperm(4);              % index pointer array
        rot = (0:1:npop-1);               % rotating index array (size NP)
        
        a1  = randperm(npop);             % shuffle locations of vectors
        rt = rem(rot+ind(1),npop);        % rotate indices by ind(1) positions
        a2  = a1(rt+1);                 % rotate vector locations
        rt = rem(rot+ind(2),npop);
        a3  = a2(rt+1);
        
        e = rand(1,lx) < CR;          % all random numbers < 0.5 are 1, 0 otherwise
        
        strategy='rand';%'best';
        switch strategy
            case 'best'
                af=randperm(nf);
                for j=1:nf
%                    [dummy,idf(j)]=min(f(:,j));
                    [dummy,idf(j)]=min(memories(:,lx+j));                                        
                end                
                %xbest=x(idf(af(1)),:)                
                xbest=memories(idf(af(1)),1:lx);                
                dx =((xbest-x(i,:)) + F*(x(a1(i),:) - x(a2(i),:)));       % differential variation
            case 'rand'
                dx =(x(a3(i),:)-x(i,:) + F*(x(a1(i),:) - x(a2(i),:)));       % differential variation
        end
        xtrial(i,:)=e.*dx+x(i,:);
        
        
        % Double check feasibility error
        for jj=1:lx
            if((xtrial(i,jj)-vlb(jj))<0||(xtrial(i,jj)-vub(jj))>0)
                xtrial(i,jj)=rand*(vub(jj)-vlb(jj))+vlb(jj);
            end
        end
        
        % evaluates move
        if norm(xtrial(i,:)-x(i,:))>0           
            if cp==0
                ftrial(i,:)=func(xtrial(i,:),arg{:});
                maxC(i)=0;
                if ~all(ftrial(i,:)>=f(i,:))
                    discarded(i).f=[discarded(i).f; ftrial(i,:)];
                    discarded(i).x=[discarded(i).x; xtrial(i,:)];
                    discarded(i).c=[discarded(i).c; maxC(i)];
                    %---------------------------------------------------------------                        
                    %  CAREFUL THIS HAS TO BE TESTED
                    if cpat
                    vtrial(i,:)=xtrial(i,:)-x(i,:);
                    end
                    %---------------------------------------------------------------                        
                    impr=1;
                end
            else
                
                [ftrial(i,:),ctrial(i,:)]=func(xtrial(i,:),arg{:});
                maxC(i)=max(ctrial(i,:));
                if cid(i)<=0&&maxC(i)>0
                    ftrial(i,:)=f(i,:)+maxC(i);
                end
                
                if (~all(ftrial(i,:)>=f(i,:))&&maxC(i)<=0)||(cid(i)>0&&maxC(i)<cid(i))
                    discarded(i).f=[discarded(i).f; ftrial(i,:)];
                    discarded(i).x=[discarded(i).x; xtrial(i,:)];
                    discarded(i).c=[discarded(i).c; maxC(i)];
                    %---------------------------------------------------------------                        
                    %  CAREFUL THIS HAS TO BE TESTED
                    if cpat
                    vtrial(i,:)=xtrial(i,:)-x(i,:);
                    end
                    %---------------------------------------------------------------                        
                    
                    impr=1;
                end
            end
            
            nfeval=nfeval+1;
            
        else
            ftrial(i,:)=f(i,:);
        end
    end
    %%
    %   PATTERN SEARCH
    if impr==0
 
         xtrial(i,:)=cpat*xtrial(i,:)+(1-cpat)*x(i,:);
%         ftrial(i,:)=f(i,:);
        xsample=[];
        fsample=[];
        
        for j=1:n_ids
            xtrial(i,ids(j))=xtrial(i,ids(j))+r(j)*rho(i,1)*Delta(ids(j));
            if xtrial(i,ids(j))<vlb(ids(j))
                xtrial(i,ids(j))=vlb(ids(j));
            elseif xtrial(i,ids(j))>vub(ids(j))
                xtrial(i,ids(j))=vub(ids(j));
            end
            
            if norm(xtrial(i,:)-x(i,:))>0
                
                if cp==0
                    ftrial(i,:)=func(xtrial(i,:),arg{:});
                    maxC(i)=0;
                    if ~all(ftrial(i,:)>=f(i,:))
                        discarded(i).f=[discarded(i).f; ftrial(i,:)];
                        discarded(i).x=[discarded(i).x; xtrial(i,:)];
                        discarded(i).c=[discarded(i).c; maxC(i)];
                        impr=1;
                    end
                else
                    [ftrial(i,:),ctrial(i,:)]=func(xtrial(i,:),arg{:});
                    maxC(i)=max(ctrial(i,:));
                    if cid(i)<=0&&maxC(i)>0
                        ftrial(i,:)=f(i,:)+maxC(i);
                    end
                    %                     xtrial(i,:)
                    %                     ftrial(i,:)
                    if (~all(ftrial(i,:)>=f(i,:))&&maxC(i)<=0)||(cid(i)>0&&maxC(i)<cid(i))
                        discarded(i).f=[discarded(i).f; ftrial(i,:)];
                        discarded(i).x=[discarded(i).x; xtrial(i,:)];
                        discarded(i).c=[discarded(i).c; maxC(i)];
                        impr=1;
                    end
                end
                nfeval=nfeval+1;
                
            end
            
            
            if (all(ftrial(i,:)>=f(i,:))&&maxC(i)<=0)||(maxC(i)>0&&maxC(i)>cid(i))
                
                xsample=[xsample;xtrial(i,:)];
                fsample=[fsample;ftrial(i,:)];
                
                if r(j)>0
                    rr=rand;
                else
                    rr=-rand;
                end
                xtrial(i,ids(j))=x(i,ids(j))-rr*rho(i,1)*Delta(ids(j));
                if xtrial(i,ids(j))<vlb(ids(j))
                    xtrial(i,ids(j))=vlb(ids(j));
                elseif xtrial(i,ids(j))>vub(ids(j))
                    xtrial(i,ids(j))=vub(ids(j));
                end
                
                if norm(xtrial(i,:)-x(i,:))>0
                    
                    if cp==0
                        ftrial(i,:)=func(xtrial(i,:),arg{:});
                        maxC(i)=0;
                        if ~all(ftrial(i,:)>=f(i,:))
                            discarded(i).f=[discarded(i).f; ftrial(i,:)];
                            discarded(i).x=[discarded(i).x; xtrial(i,:)];
                            discarded(i).c=[discarded(i).c; maxC(i)];
                            impr=1;
                        end
                    else
                        [ftrial(i,:),ctrial(i,:)]=func(xtrial(i,:),arg{:});
                        maxC(i)=max(ctrial(i,:));
                        if cid(i)<=0&&maxC(i)>0
                            ftrial(i,:)=f(i,:)+maxC(i);
                        end
                        %                         xtrial(i,:)
                        %                         ftrial(i,:)
                        if (~all(ftrial(i,:)>=f(i,:))&&maxC(i)<=0)||(cid(i)>0&&maxC(i)<cid(i))
                            discarded(i).f=[discarded(i).f; ftrial(i,:)];
                            discarded(i).x=[discarded(i).x; xtrial(i,:)];
                            discarded(i).c=[discarded(i).c; maxC(i)];
                            impr=1;
                        end
                    end
                    
                    nfeval=nfeval+1;
                    
                    if (all(ftrial(i,:)>=f(i,:))&&maxC(i)<=0)||(maxC(i)>0&&maxC(i)>cid(i))
                        
                        xsample=[xsample;xtrial(i,:)];
                        fsample=[fsample;ftrial(i,:)];
                        xtrial(i,ids(j))=x(i,ids(j));
                        ftrial(i,:)=f(i,:);
                    elseif (~all(ftrial(i,:)>=f(i,:))&&maxC(i)<=0)||(cid(i)>0&&maxC(i)<cid(i))
                        impr=1;
                        vtrial(i,:)=xtrial(i,:)-x(i,:);
                    end
                else
                    xtrial(i,ids(j))=x(i,ids(j));
                    ftrial(i,:)=f(i,:);
                end
                
            elseif (~all(ftrial(i,:)>=f(i,:))&&maxC(i)<=0)||(cid(i)>0&&maxC(i)<cid(i))
                impr=1;
                vtrial(i,:)=xtrial(i,:)-x(i,:);
            end
            if impr==1           
                break
            end
        end
        
    end
%     impr
%     pause

    if impr==0 
%         kkk=0;
%         for kk=1:2:length(fsample(:,1))-1
%             kkk=kkk+1;
%             for l=1:mfit
%                 dxsamp(kkk,:,l)=quad1D([x(i,:);xsample(kk,:);xsample(kk+1,:)],[f(i,l);fsample(kk,l);fsample(kk+1,l)],0);
%             end
%         end

%         [i pigr(act_subpr(id_pop_act_subpr==i))]
%         pause
        if  any(i==id_pop_act_subpr)&&MBHflag>0&&pigr(act_subpr(id_pop_act_subpr==i))<0.8%<- TO BE REVISED WITH A PROPPER PARAMETER TO BE SET BY THE USER
            niter=0;
            fu=1e37;
            fcurrent=g_MBHfun(x(i,:),func,lambda(act_subpr(id_pop_act_subpr==i),:),x(i,:),arg);
            while fu>=fcurrent&&niter<MBHflag                
                niter=niter+1
                pause
                foptionsNLP=optimset('Display','none','MaxFunEvals',100*lx,'LargeScale','off','Algorithm','sqp','TolFun',1e-8);
                rvlb=x(i,:)-rho(i,1)*Delta;
                rvub=x(i,:)+rho(i,1)*Delta;
                rvlb=max([rvlb;vlb]);
                rvub=min([rvub;vub]);
                xtrial(i,:)=x(i,:)+(2*rand-1)*rho(i,1)*Delta;
                xtrial(i,:)=max([xtrial(i,:);rvlb]);
                xtrial(i,:)=min([xtrial(i,:);rvub]);
                Z=x(i,:);
                [u,fu,exitflag,output]=fmincon(@(x)g_MBHfun(x,func,lambda(act_subpr(id_pop_act_subpr==i),:),Z,arg),xtrial(i,:),[],[],[],[],rvlb,rvub,@(x)c_MBHfun(x,func,lambda(act_subpr(id_pop_act_subpr==i),:),Z,arg),foptionsNLP);
                nfeval=nfeval+output.funcCount;
                if fu<fcurrent
                    xtrial(i,:)=u;
                    ftrial(i,:)=func(xtrial(i,:),arg{:});                
                    nfeval=nfeval+1;                    
                    impr=1;
                end
            end
        end
      
    
    end
    if impr==0 || all(ftrial(i,:)<f(i,:))
        rho(i,1)=rho(i,1)*contr_ratio;
 %       rho
 %       pause
        rho(i,2)=rho(i,2)+1;
        if rho(i,2)>5
            rho(i,1)=rhoini;
            rho(i,2)=0;
        end
        
     else
         
         if all(ftrial(i,:)<f(i,:))
             
             rho(i,1)=rho(i,1)/params.contr_ratio;
             rho(i,2)=rho(i,2)-1;
             
             if rho(i,2)<0                                               % unless it has contracted 3 times, in that case restart (seems to never kick in...)
                 
                 rho(i,1)=params.rhoini;
                 rho(i,2)=0;
                 
             end
             
         end
            
    end
    
    succmoves(i) = nfeval-nfbeg;
    
end

return
