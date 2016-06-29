function [xtrial,vtrial,ftrial,maxC,nfeval,discarded,rho,patdirs,MBH_positions,MADS_dirs,loc_opt]=explore2(memories,x,v,f,cid,nfeval,lambda,act_subpr,id_pop_act_subpr,z,zstar,rho,patdirs,pigr,MBH_positions,MADS_dirs,local_only,params)
%
%  [x,v,f,maxC,nfeval,discarded,rho]=explore(agents,memories,x,v,f,cid,nfeval,rho,params)
%
%
%  INPUT
%           memories                : archive if previous solutions
%           x                       : positions of agents
%           v                       : velocities of agents
%           f                       : objective values of agents
%           cid                     : max constraint violation for each agent
%           nfeval                  : number of function evaluations
%           rho                     : exploration  radius (stays as a separate input since it can change)
%
%           params : structure of parameters
%               'func'              : function handle to objective function
%               'arg'               : arguments of params.func
%               'F'                 : perturbation factor for in DE
%               'CR'                : crossover ratio for DE
%               'vlb'               : lower boundary
%               'vub'               : upper boundary
%               'coord_ratio'       : ratio of (max number of) coordinates to be explored in pattern search
%               'contr_ratio'       : contraction ratio for rho
%               'rhoini'            : initial search radius
%               'cp'                : flag for constrained or unconstrained problem
%               'ncon'              : number of constraints
%               'DE_strategy'       : strategy to use for DE
%
%  OUTPUT
%           x                       : updated position of agents
%           v                       : updated velocity of agents
%           f                       : updated function value of agents
%           maxC                    : updated max constraint violation of agents
%           nfeval                  : updated number of function evaluations
%           discarded               : structure containing successfully improved agents (why call it discarded?!)
%           rho                     : updated search radius

% CHANGE LOG
% Created by: Massimiliano Vasile 2010
% Revised, cleaned and optimized by: Lorenzo A. Ricciardi 2015


%% INITIALIZATION

xtrial=x;                                                                   % vector of position of trial agents
vtrial=0*x;                                     % vector of velocities of trial agents
ctrial=x;                                                                   % vector of contstraints of trial agents
ftrial=f;                                                                   % vector of objective function of trials agents
maxC=zeros(size(x,1),1);                                               % vector containing the max constraint violation of trial agents
discarded.f=[];                                                             % structure containing objective function values of "discarded" agents
discarded.x=[];                                                             % structure containing position of "discarded" agents
discarded.c=[];                                                             % structure containing constraint violation of "discarded" agents
Delta=(params.vub-params.vlb)/2;                                            % half span of search space
[npop,nf]=size(f);                                                          % number of agents within the population and dimensionality of criteria space (dimension of objective function f)
lx = length(params.vlb);                                                    % problem dimensionality
n_ids=ceil(sum(params.vars_to_opt)*params.coord_ratio);                                          % n_ids is equal to the number of dimensions which will actually be scanned (not automatically all of them, if params.coord_ratio <1 )
n_a = size(x,1);
mfit = size(f,2);
foptionsNLP=optimset('Display','off','MaxFunEvals',1000*length(params.vlb),'MaxIter',1000*length(params.vlb),'TolFun',1e-6,'Algorithm','sqp','MaxSQPIter', 10*length(params.vlb));
loc_opt = zeros(n_a,1);

%% MAIN LOOP
impr = zeros(n_a,1);
for i=1:n_a                                                                 % for all agents which will perform LOCAL ACTIONS
    
    %% LOOP VARIABLES RESET
    
    ids=randperm(length(params.id_vars_to_opt));                            % create a vector from 1 to n_ids, and shuffle it's elements
    ids=params.id_vars_to_opt(ids);
    r=-1+2*rand(1,n_ids);                                                   % create a vector of random numbers between -1 and 1, with n_ids elements (so as the num of coordinates which will be employed in the local search)
    
    %% INERTIA (continue along the direction of the previous improvement)
    
    if norm(v(i,:))~=0 && ~isnan(norm(v(i,:))) && ~local_only                                             % if velocity of that agent is not null (velocity is supplied as a parameter, so has been computed in previous steps). THIS IS ALWAYS FALSE DURING THE FIRST CALL OF EXPLORE
        
        alph=rand();                                                 % create another random number alpha between 0 and 1
        
        %[alph]=alpha_clip(x(i,:),alph,v(i,:),params.vlb,params.vub);                        % new, faster way of clipping_alpha (still doesn't consider alternative behaviour, can be modified inside function)
        [alph,v(i,params.id_vars_to_opt)]=alpha_clip2(x(i,params.id_vars_to_opt),alph,v(i,params.id_vars_to_opt),params.vlb(params.id_vars_to_opt),params.vub(params.id_vars_to_opt));                        % new, faster way of clipping_alpha (still doesn't consider alternative behaviour, can be modified inside function)
        
        xtrial(i,:) = x(i,:);
        
        xtrial(i,params.id_vars_to_opt)=alph*rho(i,1)*v(i,params.id_vars_to_opt)+x(i,params.id_vars_to_opt);
        xtrial(i,params.id_vars_to_opt)=xtrial(i,params.id_vars_to_opt).*(xtrial(i,params.id_vars_to_opt)>params.vlb(params.id_vars_to_opt))+params.vlb(params.id_vars_to_opt).*(xtrial(i,params.id_vars_to_opt)<=params.vlb(params.id_vars_to_opt)); % hard clipping
        xtrial(i,params.id_vars_to_opt)=xtrial(i,params.id_vars_to_opt).*(xtrial(i,params.id_vars_to_opt)<params.vub(params.id_vars_to_opt))+params.vub(params.id_vars_to_opt).*(xtrial(i,params.id_vars_to_opt)>=params.vub(params.id_vars_to_opt)); % hard clipping
        
        if norm(xtrial(i,:)-x(i,:))>0
            
            %% EVALUATE MOVE
            
            if params.cp==0                                                 % if unconstrained problem
                
                if params.optimal_control==1
                    
                    [ftrial(i,:),xcorr]=params.func(xtrial(i,:),params.arg{:});         % eval f on given position
                    xtrial(i,:) = xcorr;                                                % update initial random guess with feasible one
                    
                else
                    
                    ftrial(i,:)=params.func(xtrial(i,:),params.arg{:});         % eval f on given position
                    
                end
                
                maxC(i)=0;
                
                %if (all(ftrial(i,:)<=f(i,:))&&(norm(xtrial(i,:)-x(i,:))>0))||(any(i==id_pop_act_subpr)&&g_fun(ftrial(i,:),lambda(act_subpr(id_pop_act_subpr==i),:),z)<g_fun(f(i,:),lambda(act_subpr(id_pop_act_subpr==i),:),z)) % if agent performing local search has moved and it's objective function value has improved OR if this agent was selected to solve a sub-problem and it's partial objective function value is better than previous one
                
                if any(ftrial(i,:)<f(i,:))                          % if any component of the new trial is better than the same component of the previous x, we have an improvement!
                    
                    discarded.f=[discarded.f; ftrial(i,:)];           % put the values computed with x_trial in the discarded structure
                    discarded.x=[discarded.x; xtrial(i,:)];
                    discarded.c=[discarded.c; maxC(i)];
                    
                    if (all(ftrial(i,:)<=f(i,:))&&(norm(xtrial(i,:)-x(i,:))>0))||(any(i==id_pop_act_subpr)&&g_fun(ftrial(i,:),lambda(act_subpr(id_pop_act_subpr==i),:),z,zstar)<g_fun(f(i,:),lambda(act_subpr(id_pop_act_subpr==i),:),z,zstar)) % if agent performing local search has moved and it's objective function value has improved OR if this agent was selected to solve a sub-problem and it's partial objective function value is better than previous one
                        
                        impr(i)=1;                                                 % and say an improvement has been made, terinating the cycle
                        
                    end
                    
                    % IDEA, BETTER TO CONTINUE THE LOOP THAN CHECK THIS VAR ONWARDS, SAVES TIME AND READABILITY
                    %else
                    
                    %vtrial(i,:) = 0*vtrial(i,:);
                    
                end
                
            else                                                            % if unconstrained problem
                
                [ftrial(i,:),ctrial(i,:)]=params.func(xtrial(i,:),params.arg{:});   % eval f on given position
                maxC(i)=max(ctrial(i,:));                                           % and relative constraints violations
                
                if cid(i)<=0&&maxC(i)>0                                     % if previous position of this agent was in feasible region and now is not
                    
                    ftrial(i,:)=f(i,:)+maxC(i);                             % penalize f adding the constraint violation
                    
                end
                
                if (any(ftrial(i,:)<f(i,:))&&maxC(i)<=0)||(cid(i)>0&&maxC(i)<cid(i))	% as above: if at least 1 obective improves and respects all constraints OR at least one constraint improves from the previous iteration, keep the agent (i.e, discard it...)
                    
                    discarded.f=[discarded.f; ftrial(i,:)];
                    discarded.x=[discarded.x; xtrial(i,:)];
                    discarded.c=[discarded.c; maxC(i)];
                    
                    impr(i)=1;
                    
                end
                
            end
            
            nfeval=nfeval+1;
            
            %         else
            %             v(i,params.id_vars_to_opt)
            %             fprintf('vector clipping has screwed up things')
            %
        end
        
    end
    
    %% DDS (based on position of other agents)
    
    %     if impr==10&&any(i==id_pop_act_subpr)&&size(unique(f,'rows'),1)>2
    %
    %         %keyboard
    %         tstep=rand();
    %         [vdds,J,nu]=dds(f(1:end~=i,:),x(1:end~=i,:),maxC,f(i,:),x(i,:),lambda(act_subpr(id_pop_act_subpr==i),:),1e-3,1);
    %         dtrial=vdds;
    %
    %         %normalisation of dtrial
    %
    %         %[tstep]=alpha_clip(x(i,:),tstep,vdds,params.vlb,params.vub);                        % new, faster way of clipping_alpha (still doesn't consider alternative behaviour, can be modified inside function)
    %         [tstep,vdds]=alpha_clip2(x(i,:),tstep,vdds,params.vlb,params.vub);                        % new, faster way of clipping_alpha (still doesn't consider alternative behaviour, can be modified inside function)
    %
    %         xtrial(i,:)=tstep*vdds+x(i,:);
    %         xtrial(i,:)=xtrial(i,:).*(xtrial(i,:)>params.vlb)+params.vlb.*(xtrial(i,:)<=params.vlb)+params.vub.*(xtrial(i,:)>params.vub);  % Strong clipping, sometimes, xtrial was out of bounds by an amount close to machine eps... This led to LARGE problems (NaNs, complex numbers, etc)
    %
    %         %dtrial=vdds*tstep;
    %         %xtrial(i,:)=x(i,:)+dtrial;
    %
    %         % evaluates move
    %         if norm(xtrial(i,:)-x(i,:))>0
    %             if params.cp==0
    %                 ftrial(i,:)=params.func(xtrial(i,:),params.arg{:});
    %                 maxC(i)=0;
    %                 if ~all(ftrial(i,:)>=f(i,:))
    %                     discarded.f=[discarded.f; ftrial(i,:)];
    %                     discarded.x=[discarded.x; xtrial(i,:)];
    %                     discarded.c=[discarded.c; maxC(i)];
    %                     %---------------------------------------------------------------
    %                     %  CAREFUL THIS HAS TO BE TESTED
    %                     if params.v
    %                         %vtrial(i,:)=xtrial(i,:)-x(i,:);
    %                     end
    %                     %---------------------------------------------------------------
    %                     impr=1;
    %                 end
    %             else
    %
    %                 [ftrial(i,:),ctrial(i,:)]=params.func(xtrial(i,:),params.arg{:});
    %                 maxC(i)=max(ctrial(i,:));
    %                 if cid(i)<=0&&maxC(i)>0
    %                     ftrial(i,:)=f(i,:)+maxC(i);
    %                 end
    %
    %                 if (~all(ftrial(i,:)>=f(i,:))&&maxC(i)<=0)||(cid(i)>0&&maxC(i)<cid(i))
    %                     discarded.f=[discarded.f; ftrial(i,:)];
    %                     discarded.x=[discarded.x; xtrial(i,:)];
    %                     discarded.c=[discarded.c; maxC(i)];
    %                     %---------------------------------------------------------------
    %                     %  CAREFUL THIS HAS TO BE TESTED
    %                     if params.v
    %                         %vtrial(i,:)=xtrial(i,:)-x(i,:);
    %                     end
    %                     %---------------------------------------------------------------
    %
    %                     impr=1;
    %                 end
    %             end
    %
    %             nfeval=nfeval+1;
    %
    %         else
    %             ftrial(i,:)=f(i,:);
    %         end
    %
    %     end
    
    %% PATTERN SEARCH
    
    if impr(i)==0 && strcmp(params.pat_search_strategy,'standard') && ~local_only
        
        for j=1:n_ids                                                       % for every parameter (remember that they are shuffled every time)
            
            xtrial(i,:) = x(i,:);
            
            xtrial(i,ids(j))=xtrial(i,ids(j))+r(j)*rho(i,1)*Delta(ids(j));  % take a step in that direction, with random amplitude between 0 and current max allowed by rho (either 1/2 of total span, 1/4 or 1/8)
            
            if xtrial(i,ids(j))<params.vlb(ids(j))                                 % check that the step is in the feasible area
                
                xtrial(i,ids(j))=params.vlb(ids(j));
                
            elseif xtrial(i,ids(j))>params.vub(ids(j))
                
                xtrial(i,ids(j))=params.vub(ids(j));
                
            end
            
            if norm(xtrial(i,:)-x(i,:))>0                                   % if agent has moved from previous position
                
                %% EVALUATE MOVE
                
                if params.cp==0                                                    % if problem is unconstrained
                    
                    if params.optimal_control==1
                        
                        [ftrial(i,:),xcorr]=params.func(xtrial(i,:),params.arg{:});         % eval f on given position
                        xtrial(i,:) = xcorr;                                                % update initial random guess with feasible one
                        
                    else
                        
                        ftrial(i,:)=params.func(xtrial(i,:),params.arg{:});         % eval f on given position
                        
                    end
                    
                    maxC(i)=0;
                    
                    if any(ftrial(i,:)<f(i,:))                              % if there has been an improvement
                        
                        discarded.f=[discarded.f; ftrial(i,:)];       % add this solution to the discarded array
                        discarded.x=[discarded.x; xtrial(i,:)];
                        discarded.c=[discarded.c; maxC(i)];
                        
                        if (all(ftrial(i,:)<=f(i,:))&&(norm(xtrial(i,:)-x(i,:))>0))||(any(i==id_pop_act_subpr)&&g_fun(ftrial(i,:),lambda(act_subpr(id_pop_act_subpr==i),:),z,zstar)<g_fun(f(i,:),lambda(act_subpr(id_pop_act_subpr==i),:),z,zstar)) % if agent performing local search has moved and it's objective function value has improved OR if this agent was selected to solve a sub-problem and it's partial objective function value is better than previous one
                            
                            impr(i)=1;
                            
                        end
                        
                    end
                    
                else                                                        % if it's constrained, do as above
                    
                    [ftrial(i,:),ctrial(i,:)]=params.func(xtrial(i,:),params.arg{:});
                    maxC(i)=max(ctrial(i,:));
                    
                    if cid(i)<=0&&maxC(i)>0
                        
                        ftrial(i,:)=f(i,:)+maxC(i);
                        
                    end
                    
                    if (any(ftrial(i,:)<f(i,:))&&maxC(i)<=0)||(cid(i)>0&&maxC(i)<cid(i))
                        
                        discarded.f=[discarded.f; ftrial(i,:)];
                        discarded.x=[discarded.x; xtrial(i,:)];
                        discarded.c=[discarded.c; maxC(i)];
                        
                        impr(i)=1;
                        
                    end
                    
                end
                
                nfeval=nfeval+1;
                
            end
            
            if impr(i) == 0
                
                if r(j)>0                                                   % if previous rand number was positive
                    
                    rr=rand;                                                % generate a new positive one
                    
                else                                                        % otherwise
                    
                    rr=-rand;                                               % generate a new negative one
                    
                end
                
                xtrial(i,:) = x(i,:);
                
                xtrial(i,ids(j))=x(i,ids(j))-rr*rho(i,1)*Delta(ids(j));     % the NEW trial position is on the OPPOSITE side of the previously sampled one (notice +r before and -rr now)
                
                if xtrial(i,ids(j))<params.vlb(ids(j))                             % always check bounds
                    
                    xtrial(i,ids(j))=params.vlb(ids(j));
                    
                elseif xtrial(i,ids(j))>params.vub(ids(j))
                    
                    xtrial(i,ids(j))=params.vub(ids(j));
                    
                end
                
                if norm(xtrial(i,:)-x(i,:))>0                               % if agent has moved from initial position
                    
                    %% EVALUATE MOVE
                    
                    if params.cp==0                                                % if problem is unconstrained
                        
                        if params.optimal_control==1
                            
                            [ftrial(i,:),xcorr]=params.func(xtrial(i,:),params.arg{:});         % eval f on given position
                            xtrial(i,:) = xcorr;                                                % update initial random guess with feasible one
                            
                        else
                            
                            ftrial(i,:)=params.func(xtrial(i,:),params.arg{:});         % eval f on given position
                            
                        end
                        
                        maxC(i)=0;
                        
                        if any(ftrial(i,:)<f(i,:))                          % if it's good, save it
                            
                            discarded.f=[discarded.f; ftrial(i,:)];
                            discarded.x=[discarded.x; xtrial(i,:)];
                            discarded.c=[discarded.c; maxC(i)];
                            
                            if (all(ftrial(i,:)<=f(i,:))&&(norm(xtrial(i,:)-x(i,:))>0))||(any(i==id_pop_act_subpr)&&g_fun(ftrial(i,:),lambda(act_subpr(id_pop_act_subpr==i),:),z,zstar)<g_fun(f(i,:),lambda(act_subpr(id_pop_act_subpr==i),:),z,zstar)) % if agent performing local search has moved and it's objective function value has improved OR if this agent was selected to solve a sub-problem and it's partial objective function value is better than previous one
                                
                                impr(i)=1;
                                
                            end
                            
                        end
                        
                    else                                                    % if problem is unconstrained
                        
                        [ftrial(i,:),ctrial(i,:)]=params.func(xtrial(i,:),params.arg{:}); % evaluate new position
                        maxC(i)=max(ctrial(i,:));
                        
                        if cid(i)<=0&&maxC(i)>0
                            
                            ftrial(i,:)=f(i,:)+maxC(i);
                            
                        end
                        
                        if (any(ftrial(i,:)<f(i,:))&&maxC(i)<=0)||(cid(i)>0&&maxC(i)<cid(i))    % if it's good, keep it
                            
                            discarded.f=[discarded.f; ftrial(i,:)];
                            discarded.x=[discarded.x; xtrial(i,:)];
                            discarded.c=[discarded.c; maxC(i)];
                            
                            impr(i)=1;
                            
                        end
                        
                    end
                    
                    nfeval=nfeval+1;
                    
                end
                
            end
            
            if impr(i)==1
                
                vtrial(i,:) = xtrial(i,:)-x(i,:);
                
                break
                
            end
            
        end
        
    end
    
    %% ONE SIDED PATTERN SEARCH WITH STEP TRACKING
    
    if impr(i)==0 && strcmp(params.pat_search_strategy,'tracking') && ~local_only
        
        for j=1:n_ids                                                       % for every parameter (remember that they are shuffled every time)
            
            xtrial(i,:) = x(i,:);
            
            % reset available directions if there are no more available
            if length(patdirs(i).avail)<1
                
                patdirs(i).avail = 1:lx;
                
            end
            
            % pick available random direction
            ids = patdirs(i).avail(randperm(length(patdirs(i).avail)));
            id = ids(1);
            
            % remove selected direction from the list of available ones
            patdirs(i).avail(patdirs(i).avail==id) = [];
            
            if (xtrial(i,id)==params.vlb(id) && r(j)<0) || (xtrial(i,id)==params.vub(id) && r(j)>0)
                
                r(j) = -r(j);
                
            end
            
            xtrial(i,id)=xtrial(i,id)+r(j)*rho(i,1)*Delta(id);  % take a step in that direction, with random amplitude between 0 and current max allowed by rho (either 1/2 of total span, 1/4 or 1/8)
            
            if xtrial(i,id)<params.vlb(id)                                 % check that the step is in the feasible area
                
                xtrial(i,id)=params.vlb(id);
                
            elseif xtrial(i,id)>params.vub(id)
                
                xtrial(i,id)=params.vub(id);
                
            end
            
            if norm(xtrial(i,:)-x(i,:))>0                                   % if agent has moved from previous position
                
                %% EVALUATE MOVE
                
                if params.cp==0                                                    % if problem is unconstrained
                    
                    if params.optimal_control==1
                        
                        [ftrial(i,:),xcorr]=params.func(xtrial(i,:),params.arg{:});         % eval f on given position
                        xtrial(i,:) = xcorr;                                                % update initial random guess with feasible one
                        
                    else
                        
                        ftrial(i,:)=params.func(xtrial(i,:),params.arg{:});         % eval f on given position
                        
                    end
                    
                    maxC(i)=0;
                    
                    if any(ftrial(i,:)<f(i,:))                              % if there has been an improvement
                        
                        discarded.f=[discarded.f; ftrial(i,:)];       % add this solution to the discarded array
                        discarded.x=[discarded.x; xtrial(i,:)];
                        discarded.c=[discarded.c; maxC(i)];
                        
                        if (all(ftrial(i,:)<=f(i,:))&&(norm(xtrial(i,:)-x(i,:))>0))||(any(i==id_pop_act_subpr)&&g_fun(ftrial(i,:),lambda(act_subpr(id_pop_act_subpr==i),:),z,zstar)<g_fun(f(i,:),lambda(act_subpr(id_pop_act_subpr==i),:),z,zstar)) % if agent performing local search has moved and it's objective function value has improved OR if this agent was selected to solve a sub-problem and it's partial objective function value is better than previous one
                            
                            impr(i)=1;
                            
                        end
                        
                    end
                    
                else                                                        % if it's constrained, do as above
                    
                    [ftrial(i,:),ctrial(i,:)]=params.func(xtrial(i,:),params.arg{:});
                    maxC(i)=max(ctrial(i,:));
                    
                    if cid(i)<=0&&maxC(i)>0
                        
                        ftrial(i,:)=f(i,:)+maxC(i);
                        
                    end
                    
                    if (any(ftrial(i,:)<f(i,:))&&maxC(i)<=0)||(cid(i)>0&&maxC(i)<cid(i))
                        
                        discarded.f=[discarded.f; ftrial(i,:)];
                        discarded.x=[discarded.x; xtrial(i,:)];
                        discarded.c=[discarded.c; maxC(i)];
                        
                        impr(i)=1;
                        
                    end
                    
                end
                
                nfeval=nfeval+1;
                
            end
            
            if impr(i)==1
                
                vtrial(i,:) = xtrial(i,:)-x(i,:);
                
                break
                
            end
            
        end
        
    end
    
    %% RANDOM DIRECTION PATTERN SEARCH
    
    if impr(i)==0 && strcmp(params.pat_search_strategy,'random') && ~local_only
        
        for j=1:n_ids                                                       % for every parameter (remember that they are shuffled every time)
            
            xtrial(i,:) = x(i,:);
            
            % random perturbation of size rho(i,1)
            pert = rand(1,length(params.id_vars_to_opt));
            alph = rand()*rho(i,1);
            
            [alph,pert] = alpha_clip2(x(i,params.id_vars_to_opt),alph,pert,params.vlb(params.id_vars_to_opt),params.vub(params.id_vars_to_opt));
            
            xtrial(i,:) = x(i,:);
            xtrial(i,params.id_vars_to_opt)=alph*pert+x(i,params.id_vars_to_opt);
            xtrial(i,params.id_vars_to_opt)=xtrial(i,params.id_vars_to_opt).*(xtrial(i,params.id_vars_to_opt)>params.vlb(params.id_vars_to_opt))+params.vlb(params.id_vars_to_opt).*(xtrial(i,params.id_vars_to_opt)<=params.vlb(params.id_vars_to_opt)); % hard clipping
            xtrial(i,params.id_vars_to_opt)=xtrial(i,params.id_vars_to_opt).*(xtrial(i,params.id_vars_to_opt)<params.vub(params.id_vars_to_opt))+params.vub(params.id_vars_to_opt).*(xtrial(i,params.id_vars_to_opt)>=params.vub(params.id_vars_to_opt)); % hard clipping
            
            if norm(xtrial(i,:)-x(i,:))>0                                   % if agent has moved from previous position
                
                %% EVALUATE MOVE
                
                if params.cp==0                                                    % if problem is unconstrained
                    
                    if params.optimal_control==1
                        
                        [ftrial(i,:),xcorr]=params.func(xtrial(i,:),params.arg{:});         % eval f on given position
                        xtrial(i,:) = xcorr;                                                % update initial random guess with feasible one
                        
                    else
                        
                        ftrial(i,:)=params.func(xtrial(i,:),params.arg{:});         % eval f on given position
                        
                    end
                    
                    maxC(i)=0;
                    
                    if any(ftrial(i,:)<f(i,:))                              % if there has been an improvement
                        
                        discarded.f=[discarded.f; ftrial(i,:)];       % add this solution to the discarded array
                        discarded.x=[discarded.x; xtrial(i,:)];
                        discarded.c=[discarded.c; maxC(i)];
                        
                        if (all(ftrial(i,:)<=f(i,:))&&(norm(xtrial(i,:)-x(i,:))>0))||(any(i==id_pop_act_subpr)&&g_fun(ftrial(i,:),lambda(act_subpr(id_pop_act_subpr==i),:),z,zstar)<g_fun(f(i,:),lambda(act_subpr(id_pop_act_subpr==i),:),z,zstar)) % if agent performing local search has moved and it's objective function value has improved OR if this agent was selected to solve a sub-problem and it's partial objective function value is better than previous one
                            
                            impr(i)=1;
                            
                        end
                        
                    end
                    
                else                                                        % if it's constrained, do as above
                    
                    [ftrial(i,:),ctrial(i,:)]=params.func(xtrial(i,:),params.arg{:});
                    maxC(i)=max(ctrial(i,:));
                    
                    if cid(i)<=0&&maxC(i)>0
                        
                        ftrial(i,:)=f(i,:)+maxC(i);
                        
                    end
                    
                    if (any(ftrial(i,:)<f(i,:))&&maxC(i)<=0)||(cid(i)>0&&maxC(i)<cid(i))
                        
                        discarded.f=[discarded.f; ftrial(i,:)];
                        discarded.x=[discarded.x; xtrial(i,:)];
                        discarded.c=[discarded.c; maxC(i)];
                        
                        impr(i)=1;
                        
                    end
                    
                end
                
                nfeval=nfeval+1;
                
            end
            
            if impr(i)==1
                
                vtrial(i,:) = xtrial(i,:)-x(i,:);
                
                break
                
            end
            
        end
        
    end
    
    %% MADS
    
    if impr(i)==0 && strcmp(params.pat_search_strategy,'MADS') && ~local_only
        
        [MADS_dirs,D] = MADS_generate_spanning_set (MADS_dirs,rho(i,2)+1,length(params.id_vars_to_opt),rho(i,1),1);
        
        for j = 1:n_ids%size(D,2)
            
            pert = D(:,j)';
            alph = rho(i,1);
            
            [alph,pert] = alpha_clip2(x(i,params.id_vars_to_opt),alph,pert,params.vlb(params.id_vars_to_opt),params.vub(params.id_vars_to_opt));
            
            xtrial(i,:) = x(i,:);
            
            xtrial(i,params.id_vars_to_opt)=alph*pert+x(i,params.id_vars_to_opt);
            xtrial(i,params.id_vars_to_opt)=xtrial(i,params.id_vars_to_opt).*(xtrial(i,params.id_vars_to_opt)>params.vlb(params.id_vars_to_opt))+params.vlb(params.id_vars_to_opt).*(xtrial(i,params.id_vars_to_opt)<=params.vlb(params.id_vars_to_opt)); % hard clipping
            xtrial(i,params.id_vars_to_opt)=xtrial(i,params.id_vars_to_opt).*(xtrial(i,params.id_vars_to_opt)<params.vub(params.id_vars_to_opt))+params.vub(params.id_vars_to_opt).*(xtrial(i,params.id_vars_to_opt)>=params.vub(params.id_vars_to_opt)); % hard clipping
            
            if norm(xtrial(i,:)-x(i,:))>0                                   % if agent has moved from previous position
                
                %% EVALUATE MOVE
                
                if params.cp==0                                                    % if problem is unconstrained
                    
                    if params.optimal_control==1
                        
                        [ftrial(i,:),xcorr]=params.func(xtrial(i,:),params.arg{:});         % eval f on given position
                        xtrial(i,:) = xcorr;                                                % update initial random guess with feasible one
                        
                    else
                        
                        ftrial(i,:)=params.func(xtrial(i,:),params.arg{:});         % eval f on given position
                        
                    end
                    
                    maxC(i)=0;
                    
                    if any(ftrial(i,:)<f(i,:))                              % if there has been an improvement
                        
                        discarded.f=[discarded.f; ftrial(i,:)];       % add this solution to the discarded array
                        discarded.x=[discarded.x; xtrial(i,:)];
                        discarded.c=[discarded.c; maxC(i)];
                        
                        if (all(ftrial(i,:)<=f(i,:))&&(norm(xtrial(i,:)-x(i,:))>0))||(any(i==id_pop_act_subpr)&&g_fun(ftrial(i,:),lambda(act_subpr(id_pop_act_subpr==i),:),z,zstar)<g_fun(f(i,:),lambda(act_subpr(id_pop_act_subpr==i),:),z,zstar)) % if agent performing local search has moved and it's objective function value has improved OR if this agent was selected to solve a sub-problem and it's partial objective function value is better than previous one
                            
                            impr(i)=1;
                            vtrial(i,:) = xtrial(i,:)-x(i,:);
                            
                        end
                        
                    end
                    
                else                                                        % if it's constrained, do as above
                    
                    [ftrial(i,:),ctrial(i,:)]=params.func(xtrial(i,:),params.arg{:});
                    maxC(i)=max(ctrial(i,:));
                    
                    if cid(i)<=0&&maxC(i)>0
                        
                        ftrial(i,:)=f(i,:)+maxC(i);
                        
                    end
                    
                    if (any(ftrial(i,:)<f(i,:))&&maxC(i)<=0)||(cid(i)>0&&maxC(i)<cid(i))
                        
                        discarded.f=[discarded.f; ftrial(i,:)];
                        discarded.x=[discarded.x; xtrial(i,:)];
                        discarded.c=[discarded.c; maxC(i)];
                        
                        impr(i)=1;
                        
                    end
                    
                end
                
                nfeval=nfeval+1;
                
            end
            
            if impr(i)==1
                
                vtrial(i,:) = xtrial(i,:)-x(i,:);
                
                break
                
            end
            
        end
        
    end
    
    %% DIFFERENTIAL (take the direction suggested by other agents)
    
    if impr(i)==0 && ~local_only
        
        %redo =1;
        %cnt = 0;
        
        %while redo %&& cnt<100
        
        %    cnt = cnt+1;
        %    if cnt>100
        %        keyboard
        %    end
        qq = randperm(npop);
        %    qq(qq==i)=[];
        a1 = qq(1);
        a2 = qq(2);
        a3 = qq(3);
        
        %         % The following prototype was to perform DE either on the current
        %         % agents' position or on the archive, in a unified and equally
        %         % probable way. Not tested.
        %
        %         tmpmm = memories(:,1:lx);%x(1:end~=i,:)];
        %         %tmpmm = unique(tmpmm,'rows');
        %         %[domA,~] = dominant([memories(:,lx+1:lx+mfit);f(1:end~=i,:)],f(i,:));
        %         %tmpmm = tmpmm(domA==0,:);
        %
        %         qq = randperm(size(tmpmm,1));
        %
        %         if length(qq)>2
        %
        %             x1 = tmpmm(qq(1));
        %             x2 = tmpmm(qq(2));
        %             x3 = tmpmm(qq(3));
        %
        %         else
        %
        %             if length(qq)==2
        %
        %                 x3 = x(i);
        %                 x1 = tmpmm(qq(1));
        %                 x2 = tmpmm(qq(2));
        %
        %             else
        %
        %                 x3 = tmpmm(qq(1));
        %                 x1 = 0*x3;
        %                 x2 = 0*x3;
        %
        %             end
        %
        %         end
        
        e = rand(1,length(params.id_vars_to_opt)) < params.CR;                                                % all random numbers < 0.5 are 1, 0 otherwise
        
        switch params.DE_strategy
            
            case 'best'
                
                af=randperm(nf);                                            % af is a vector containing spanning the number of objectives, shuffled
                idf = zeros(1,nf);
                
                % detecting best agents for each objective
                for j=1:nf                                                  % for each objective
                    
                    [~,idf(j)]=min(memories(:,lx+j));                       % pick the ID of the agent that solves it best
                    
                end
                
                xbest=memories(idf(af(1)),1:lx);                            % this way, xbest is the agent which best solves an OBJECTIVE PICKED AT RANDOM
                dx =rand()*((xbest-x(i,:)) + rand()*params.F*(x(a1,:) - x(a2,:)));         % differential variation
                
            case 'rand'
                
                dx =  (x(a3,params.id_vars_to_opt)-x(i,params.id_vars_to_opt)) + params.F*(x(a1,params.id_vars_to_opt) - x(a2,params.id_vars_to_opt));
                %dx = (x3-x(i,:)) + params.F*(x1 - x2);
        end
        
        % Double check feasibility error
        
        % alternate formulation, probably makes a bit more sense... needs
        % testing
        vv = e.*dx;
        alph =rand();
        [alph,vv] = alpha_clip2(x(i,params.id_vars_to_opt),alph,vv,params.vlb(params.id_vars_to_opt),params.vub(params.id_vars_to_opt));
        
        xtrial(i,:) = x(i,:);
        xtrial(i,params.id_vars_to_opt)=alph*vv+x(i,params.id_vars_to_opt);
        xtrial(i,params.id_vars_to_opt)=xtrial(i,params.id_vars_to_opt).*(xtrial(i,params.id_vars_to_opt)>params.vlb(params.id_vars_to_opt))+params.vlb(params.id_vars_to_opt).*(xtrial(i,params.id_vars_to_opt)<=params.vlb(params.id_vars_to_opt)); % hard clipping
        xtrial(i,params.id_vars_to_opt)=xtrial(i,params.id_vars_to_opt).*(xtrial(i,params.id_vars_to_opt)<params.vub(params.id_vars_to_opt))+params.vub(params.id_vars_to_opt).*(xtrial(i,params.id_vars_to_opt)>=params.vub(params.id_vars_to_opt)); % hard clipping
        
        
        % hard clipping formulation
        
        %xtrial(i,:)=e.*dx+x(i,:);
        %xtrial(i,:)=xtrial(i,:).*(xtrial(i,:)>params.vlb)+params.vlb.*(xtrial(i,:)<=params.vlb); % hard clipping
        %xtrial(i,:)=xtrial(i,:).*(xtrial(i,:)<params.vub)+params.vub.*(xtrial(i,:)>=params.vub); % hard clipping
        
        % native forulation
        
        %        xtrial(i,:)=e.*dx+x(i,:);
        
        %        mask = (xtrial(i,:)<params.vlb)+(xtrial(i,:)>params.vub);           % mask vector containing components out of bounds
        %        xnew = mask.*rand(1,lx).*(params.vub-params.vlb)+params.vlb;        % vector of new components, where xtrial is out of bounds
        %        xtrial(i,:) = xtrial(i,:).*~mask+xnew.*mask;                        % new vector, completely in bounds
        %if any(xtrial(i,:)~=x(i,:))
        
        %    redo=0;
        
        %end
        
        %end
        % evaluates move
        
        if norm(xtrial(i,:)-x(i,:))>0
            
            %% EVALUATE MOVE
            
            if params.cp==0
                
                if params.optimal_control==1
                    
                    [ftrial(i,:),xcorr]=params.func(xtrial(i,:),params.arg{:});         % eval f on given position
                    xtrial(i,:) = xcorr;                                                % update initial random guess with feasible one
                    
                else
                    
                    ftrial(i,:)=params.func(xtrial(i,:),params.arg{:});         % eval f on given position
                    
                end
                
                
                maxC(i)=0;
                
                %if (all(ftrial(i,:)<=f(i,:))&&(norm(xtrial(i,:)-x(i,:))>0))||(any(i==id_pop_act_subpr)&&g_fun(ftrial(i,:),lambda(act_subpr(id_pop_act_subpr==i),:),z)<g_fun(f(i,:),lambda(act_subpr(id_pop_act_subpr==i),:),z)) % if agent performing local search has moved and it's objective function value has improved OR if this agent was selected to solve a sub-problem and it's partial objective function value is better than previous one
                if any(ftrial(i,:)<f(i,:))                                  % if any component of the new trial is better than the same component of the previous x, we have an improvement!
                    
                    discarded.f=[discarded.f; ftrial(i,:)];
                    discarded.x=[discarded.x; xtrial(i,:)];
                    discarded.c=[discarded.c; maxC(i)];
                    
                    if (all(ftrial(i,:)<=f(i,:))&&(norm(xtrial(i,:)-x(i,:))>0))||(any(i==id_pop_act_subpr)&&g_fun(ftrial(i,:),lambda(act_subpr(id_pop_act_subpr==i),:),z,zstar)<g_fun(f(i,:),lambda(act_subpr(id_pop_act_subpr==i),:),z,zstar)) % if agent performing local search has moved and it's objective function value has improved OR if this agent was selected to solve a sub-problem and it's partial objective function value is better than previous one
                        
                        impr(i)=1;
                        
                        if params.v==1
                            
                            vtrial(i,:) = xtrial(i,:)-x(i,:);
                            
                        end
                        
                    end
                    
                end
                
            else
                
                [ftrial(i,:),ctrial(i,:)]=params.func(xtrial(i,:),params.arg{:});
                maxC(i)=max(ctrial(i,:));
                
                if cid(i)<=0&&maxC(i)>0
                    
                    ftrial(i,:)=f(i,:)+maxC(i);
                    
                end
                
                if (any(ftrial(i,:)<f(i,:))&&maxC(i)<=0)||(cid(i)>0&&maxC(i)<cid(i))	% this reads: if at least 1 obective improves and respects all constraints OR at least one constraint improves from the previous iteration, keep the agent (i.e, discard it...)
                    
                    discarded.f=[discarded.f; ftrial(i,:)];
                    discarded.x=[discarded.x; xtrial(i,:)];
                    discarded.c=[discarded.c; maxC(i)];
                    
                    impr(i)=1;
                    
                end
                
            end
            
            nfeval=nfeval+1;
            
        else
            
            ftrial(i,:)=f(i,:);
            
        end
        
    end
    
    %% DDS step
    
    %     if impr==100&&any(i==id_pop_act_subpr)
    %
    %         [vdds,J,nu]=dds(fsample,xsample,discarded.c(i),f0trial(i,:),x0trial(i,:),lambda(act_subpr(id_pop_act_subpr==i),:),1e-3,2);
    %         tstep=rho(i,1);
    %         delta_angle=0;
    %         tol=1e-1;
    %         step_trial=0;
    %         while delta_angle<(1-tol)&&step_trial<10&&impr==0
    %
    %             step_trial=step_trial+1;
    %             dtrial=vdds;
    %             for jj=1:lx
    %                 if((dtrial(jj)+x0trial(i,jj))>params.vub(jj))&&(dtrial(jj)~=0)
    %                     tstep=min([tstep,abs((params.vub(jj)-x0trial(i,jj))/dtrial(jj))]);
    %                 elseif ((dtrial(jj)+x0trial(i,jj))<params.vlb(jj))&&(dtrial(jj)~=0)
    %                     tstep=min([tstep,abs((x0trial(i,jj)-params.vlb(jj))/dtrial(jj))]);
    %                 end
    %             end
    %             dtrial=vdds*tstep;
    %             xtrial(i,:)=x0trial(i,:)+dtrial;
    %
    %             if any(xtrial(i,:)<params.vlb)
    %                 pause
    %             end
    %             if any(xtrial(i,:)>params.vub)
    %                 pause
    %             end
    %             if params.cp==0
    %                 ftrial(i,:)=params.func(xtrial(i,:),params.arg{:});
    %                 maxC(i)=0;
    %                 if ~all(ftrial(i,:)>=f(i,:))
    %                     discarded.f(i)=[discarded.f(i); ftrial(i,:)];
    %                     discarded.x(i)=[discarded.x(i); xtrial(i,:)];
    %                     discarded.c(i)=[discarded.c(i); maxC(i)];
    %                     impr=1;
    %
    %                     %vtrial(i,:)=xtrial(i,:)-x(i,:);
    %
    %                 end
    %             else
    %                 [ftrial(i,:),ctrial(i,:)]=params.func(xtrial(i,:),params.arg{:});
    %                 maxC(i)=max(ctrial(i,:));
    %                 if cid(i)<=0&&maxC(i)>0
    %                     ftrial(i,:)=f(i,:)+maxC(i);
    %                 end
    %                 if (~all(ftrial(i,:)>=f(i,:))&&maxC(i)<=0)||(cid(i)>0&&maxC(i)<cid(i))
    %                     discarded.f(i)=[discarded.f(i); ftrial(i,:)];
    %                     discarded.x(i)=[discarded.x(i); xtrial(i,:)];
    %                     discarded.c(i)=[discarded.c(i); maxC(i)];
    %                     impr=1;
    %
    %                     %vtrial(i,:)=xtrial(i,:)-x(i,:);
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
    
    if impr(i)==0
        
        if params.optimal_control==1
            
            % NOTE: FOR SINGLE OBJECTIVES, SMOTHED TCHEBYSHEV SCALARISATION
            % IS NOT THE BEST APPROACH, AS WE DON'T HAVE A CLUE ABOUT THE
            % UTOPIA POINT (IT IS DEFINED BY THE SINGLE OBJECTIVES!!!)
            % EITHER USE TWO SEPARATE APPROACHES, OR SETTLE DOWN FOR SUB
            % OPTIMAL SOLUTIONS (IF THE UTOPIA POINT PREDICTED IS NOT AS
            % GOOD AS THE REAL ONE, SMOOTHED TCHEBYSHEV SCALARISATION WILL
            % STOP BEFORE)
            
            if any(i==id_pop_act_subpr)&& params.MBHflag>0 && local_only%any(i==id_pop_act_subpr)&& params.MBHflag>0 && (rho(i,2)==params.max_rho_contr || local_only) %&&pigr(act_subpr(id_pop_act_subpr==i))<0.8%<- TO BE REVISED WITH A PROPPER PARAMETER TO BE SET BY THE USER
                
                if isempty(MBH_positions)
                                        
                    fprintf('Running MBH on agent %d, lambda = (',i);
                    
                    for j =1:length(lambda(act_subpr(id_pop_act_subpr==i),:))-1
                        
                        fprintf('%f ',lambda(act_subpr(id_pop_act_subpr==i),j));
                        
                    end
                    
                    fprintf('%f)\n',lambda(act_subpr(id_pop_act_subpr==i),end));
                    niter=0;
                    fu=Inf;
                    
                    fcurrent=1;
                    
                    while fu>=fcurrent&&niter<params.MBHflag
                        
                        niter=niter+1;
                        fprintf('Iter %d\n',niter)
                        rvlb=x(i,:)-rho(i,1)*Delta;
                        rvub=x(i,:)+rho(i,1)*Delta;
                        rvlb=max([rvlb;params.vlb]);
                        rvub=min([rvub;params.vub]);
                        
                        if niter>1
                            
                            xtrial(i,:)=x(i,:)+(2*rand-1)*rho(i,1)*Delta;
                            xtrial(i,:)=max([xtrial(i,:);rvlb]);
                            xtrial(i,:)=min([xtrial(i,:);rvub]);
                            [ftrial(i,:),xtrial(i,:)] = params.func(xtrial(i,:),params.arg{:});
                            
                        else
                            
                            xtrial(i,:) = x(i,:);
                            ftrial(i,:) = f(i,:);
                            
                        end
                        
                        zz = 2*z-zstar;
                        zzstar = ftrial(i,:);                        
                        
                        if act_subpr(id_pop_act_subpr==i)<=mfit
                            
                            % In this case we follow the original
                            % objective, so the utopia point is just
                            % pushed as far as possible but along the
                            % same direction
                            
                            ll = lambda(act_subpr(id_pop_act_subpr==i),:);
                            %zz = 2*z-zstar;
                            %zzstar = ftrial(i,:);
                            
                        else
                            
                            % in this case the scalarisation is looking
                            % for strongly dominant solutions, moving
                            % along the direction (1,1,1,1...,1) in
                            % criteria space. For this reason, we must
                            % move the utopia point, keeping in mind
                            % that the problem will be scaled in the
                            % nonlinear solver. Thus at first the
                            % search direction will be scaled
                            % acccording to (zstar-z), to find the
                            % "correct" utopia point, then the
                            % (1,1,1,..1) direction will be used in the
                            % gradient optimiser
                            
                            %ll = ones(1,mfit)./norm(ones(1,mfit)).*(zstar-z);
                            %ll = ll/norm(ll);
                            %ztemp = 2*z-zstar;
                            %zz = ftrial(i,:)-ll*(ftrial(i,1)-ztemp(1))/ll(1);
                            %zzstar = ftrial(i,:);
                            ll = ones(1,mfit)./norm(ones(1,mfit));
                            
                        end
                        
                        xt=[norm(ll) xtrial(i,:)];
                        [u,fu,~,output]=fmincon(@(xt) xt(1),xt,[],[],[],[],[0 params.vlb],[norm(ll) params.vub],@(xt) params.oc.smooth_scal_constr_fun(xt,ll,zz,zzstar,params),foptionsNLP);
                        nfeval=nfeval+output.funcCount+1;
                        loc_opt(i)=1;
                        
                        if params.optimal_control==0
                            
                            xtrial(i,:) = u;
                            ftrial(i,:) = params.func(u,params.arg{:});
                            
                        else
                            
                            [ftrial(i,:),xtrial(i,:)] = params.func(u(2:end),params.arg{:});
                            
                        end
                        
                        fprintf('Old solution = ');
                        
                        for j=1:mfit-1
                            
                            fprintf('%f ',f(i,j)) ;
                            
                        end
                        
                        fprintf('%f)\n',f(i,end));
                        fprintf('New solution = ');
                        
                        for j=1:mfit-1
                            
                            fprintf('%f ',ftrial(i,j)) ;
                            
                        end
                        
                        fprintf('%f)\n',ftrial(i,end));
                        
                        if any(ftrial(i,:)<f(i,:))                                  % if any component of the new trial is better than the same component of the previous x, we have an improvement!
                            
                            discarded.f=[discarded.f; ftrial(i,:)];
                            discarded.x=[discarded.x; xtrial(i,:)];
                            discarded.c=[discarded.c; maxC(i)];
                            vtrial(i,:) = xtrial(i,:)-x(i,:);
                            
                            if (all(ftrial(i,:)<=f(i,:))&&(norm(xtrial(i,:)-x(i,:))>0))||(any(i==id_pop_act_subpr)&&g_fun(ftrial(i,:),lambda(act_subpr(id_pop_act_subpr==i),:),z,zstar)<g_fun(f(i,:),lambda(act_subpr(id_pop_act_subpr==i),:),z,zstar)) % if agent performing local search has moved and it's objective function value has improved OR if this agent was selected to solve a sub-problem and it's partial objective function value is better than previous one
                                
                                impr(i)=1;
                                %MBH_positions = [MBH_positions; xtrial(i,:)];                    %avoids cascade of fmincon usage, should be done selectively if the improvement was too little for the effort. How do I measure it?
                                
                                fprintf('Success!!!\n')
                                
                                % delays next MBH as much as possible
                                rho(i,1)=params.rhoini;
                                rho(i,2)=0;
                                
                            else
                                
                                fprintf('There is a mismatch between Tchebychev scalarisation and the smooth version!\n');
                                fprintf('g_fun of old solution: %f\n',g_fun(f(i,:),lambda(act_subpr(id_pop_act_subpr==i),:),z,zstar));
                                fprintf('g_fun of trial solution: %f\n',g_fun(ftrial(i,:),lambda(act_subpr(id_pop_act_subpr==i),:),z,zstar));
                                
                            end
                            
                        else
                            
                            fprintf('Fail...\n')
                            
                        end
                        
                    end
                    
                else
                    
                    % look if this position has not already been the starting point
                    % for fmincon
                    
                    if ~any(all(repmat(x(i,:),size(MBH_positions,1),1)==MBH_positions,2))
                        
                        MBH_positions = [MBH_positions; x(i,:)];                    %add this position in the list, to avoid repeating it
                        
                        fprintf('Running MBH on agent %d, lambda = (',i);
                        
                        for j =1:length(lambda(act_subpr(id_pop_act_subpr==i),:))-1
                            
                            fprintf('%f ',lambda(act_subpr(id_pop_act_subpr==i),j));
                            
                        end
                        
                        fprintf('%f)\n',lambda(act_subpr(id_pop_act_subpr==i),end));
                        niter=0;
                        fu=Inf;
                        
                        fcurrent=1;
                        
                        while fu>=fcurrent&&niter<params.MBHflag
                            
                            niter=niter+1;
                            fprintf('Iter %d\n',niter)
                            rvlb=x(i,:)-rho(i,1)*Delta;
                            rvub=x(i,:)+rho(i,1)*Delta;
                            rvlb=max([rvlb;params.vlb]);
                            rvub=min([rvub;params.vub]);
                            
                            if niter>1
                                
                                xtrial(i,:)=x(i,:)+(2*rand-1)*rho(i,1)*Delta;
                                xtrial(i,:)=max([xtrial(i,:);rvlb]);
                                xtrial(i,:)=min([xtrial(i,:);rvub]);
                                
                            else
                                
                                xtrial(i,:) = x(i,:);
                                
                            end
                            
                            ftrial(i,:)=params.func(xtrial(i,:),params.arg{:});
                            
                            if act_subpr(id_pop_act_subpr==i)<=mfit
                                
                                % In this case we follow the original
                                % objective, so the utopia point is just
                                % pushed as far as possible but along the
                                % same direction
                                
                                ll = lambda(act_subpr(id_pop_act_subpr==i),:);
                                zz = 2*z-zstar;
                                zzstar = ftrial(i,:);
                                
                            else
                                
                                % in this case the scalarisation is looking
                                % for strongly dominant solutions, moving
                                % along the direction (1,1,1,1...,1) in
                                % criteria space. For this reason, we must
                                % move the utopia point, keeping in mind
                                % that the problem will be scaled in the
                                % nonlinear solver. Thus at first the
                                % search direction will be scaled
                                % acccording to (zstar-z), to find the
                                % "correct" utopia point, then the
                                % (1,1,1,..1) direction will be used in the
                                % gradient optimiser
                                
                                ll = ones(1,mfit)./norm(ones(1,mfit)).*(zstar-z);
                                ll = ll/norm(ll);
                                ztemp = 2*z-zstar;
                                zz = ftrial(i,:)-ll*(ftrial(i,1)-ztemp(1))/ll(1);
                                zzstar = ftrial(i,:);
                                ll = ones(1,mfit)./norm(ones(1,mfit));
                                
                            end
                            
                            xt=[norm(ll) xtrial(i,:)];
                            [u,fu,~,output]=fmincon(@(xt) xt(1),xt,[],[],[],[],[0 params.vlb],[norm(ll) params.vub],@(xt) params.oc.smooth_scal_constr_fun(xt,ll,zz,zzstar,params),foptionsNLP);
                            nfeval=nfeval+output.funcCount+1;
                            loc_opt(i)=1;
                            
                            if params.optimal_control==0
                                
                                xtrial(i,:) = u;
                                ftrial(i,:) = params.func(u,params.arg{:});
                                
                            else
                                
                                [ftrial(i,:),xtrial(i,:)] = params.func(u(2:end),params.arg{:});
                                
                            end
                            
                            fprintf('Old solution = ');
                            
                            for j=1:mfit-1
                                
                                fprintf('%f ',f(i,j)) ;
                                
                            end
                            
                            fprintf('%f)\n',f(i,end));
                            fprintf('New solution = ');
                            
                            for j=1:mfit-1
                                
                                fprintf('%f ',ftrial(i,j)) ;
                                
                            end
                            
                            fprintf('%f)\n',ftrial(i,end));
                            
                            if any(ftrial(i,:)<f(i,:))                                  % if any component of the new trial is better than the same component of the previous x, we have an improvement!
                                
                                discarded.f=[discarded.f; ftrial(i,:)];
                                discarded.x=[discarded.x; xtrial(i,:)];
                                discarded.c=[discarded.c; maxC(i)];
                                vtrial(i,:) = xtrial(i,:)-x(i,:);
                                
                                if (all(ftrial(i,:)<=f(i,:))&&(norm(xtrial(i,:)-x(i,:))>0))||(any(i==id_pop_act_subpr)&&g_fun(ftrial(i,:),lambda(act_subpr(id_pop_act_subpr==i),:),z,zstar)<g_fun(f(i,:),lambda(act_subpr(id_pop_act_subpr==i),:),z,zstar)) % if agent performing local search has moved and it's objective function value has improved OR if this agent was selected to solve a sub-problem and it's partial objective function value is better than previous one
                                    
                                    impr(i)=1;
                                    %MBH_positions = [MBH_positions; xtrial(i,:)];                    %avoids cascade of fmincon usage, should be done selectively if the improvement was too little for the effort. How do I measure it?
                                    
                                    fprintf('Success!!!\n')
                                    
                                    % delays next MBH as much as possible
                                    rho(i,1)=params.rhoini;
                                    rho(i,2)=0;
                                    
                                else
                                    
                                    fprintf('There is a mismatch between Tchebychev scalarisation and the smooth version!\n');
                                    fprintf('g_fun of old solution: %f\n',g_fun(f(i,:),lambda(act_subpr(id_pop_act_subpr==i),:),z,zstar));
                                    fprintf('g_fun of trial solution: %f\n',g_fun(ftrial(i,:),lambda(act_subpr(id_pop_act_subpr==i),:),z,zstar));
                                    
                                    
                                end
                                
                            else
                                
                                fprintf('Fail...\n')
                                
                            end
                            
                        end
                        
                    end
                    
                end
                
            end
            
        else
            
            if  any(i==id_pop_act_subpr)&&params.MBHflag>0 && rho(i,2)==params.max_rho_contr %&&pigr(act_subpr(id_pop_act_subpr==i))<0.8%<- TO BE REVISED WITH A PROPPER PARAMETER TO BE SET BY THE USER
                
                if isempty(MBH_positions)
                    
                    MBH_positions = [MBH_positions; x(i,:)];                    %add this position in the list, to avoid repeating it
                    
                    fprintf('Running MBH on agent %d, lambda = (',i);
                    
                    for j =1:length(lambda(act_subpr(id_pop_act_subpr==i),:))-1
                        
                        fprintf('%f ',lambda(act_subpr(id_pop_act_subpr==i),j));
                        
                    end
                    
                    fprintf('%f)\n',lambda(act_subpr(id_pop_act_subpr==i),end));
                    niter=0;
                    fu=Inf;
                    
                    fcurrent=g_MBHfun2(x(i,:),params.func,lambda(act_subpr(id_pop_act_subpr==i),:),z,zstar,params.arg);
                    
                    while fu>=fcurrent&&niter<params.MBHflag
                        
                        niter=niter+1;
                        fprintf('Iter %d\n',niter)
                        rvlb=x(i,:)-rho(i,1)*Delta;
                        rvub=x(i,:)+rho(i,1)*Delta;
                        rvlb=max([rvlb;params.vlb]);
                        rvub=min([rvub;params.vub]);
                        
                        if niter>1
                            
                            xtrial(i,:)=x(i,:)+(2*rand-1)*rho(i,1)*Delta;
                            xtrial(i,:)=max([xtrial(i,:);rvlb]);
                            xtrial(i,:)=min([xtrial(i,:);rvub]);
                            
                        else
                            
                            xtrial(i,:) = x(i,:);
                            
                        end
                        
                        ftrial(i,:)=params.func(xtrial(i,:),params.arg{:});
                        
                        xt=xtrial(i,:);
                        
                        [u,fu,~,output]=fmincon(@(xt)g_MBHfun2(xt,params.func,lambda(act_subpr(id_pop_act_subpr==i),:),z,zstar,params.arg),xt,[],[],[],[],params.vlb,params.vub,[],foptionsNLP);
                        nfeval=nfeval+output.funcCount+1;
                        
                        if params.optimal_control==0
                            
                            xtrial(i,:) = u;
                            ftrial(i,:) = params.func(u,params.arg{:});
                            
                        else
                            
                            [ftrial(i,:),xtrial(i,:)] = params.func(u,params.arg{:});
                            
                        end
                        
                        if any(ftrial(i,:)<f(i,:))                                  % if any component of the new trial is better than the same component of the previous x, we have an improvement!
                            
                            discarded.f=[discarded.f; ftrial(i,:)];
                            discarded.x=[discarded.x; xtrial(i,:)];
                            discarded.c=[discarded.c; maxC(i)];
                            vtrial(i,:) = xtrial(i,:)-x(i,:);
                            
                            if (all(ftrial(i,:)<=f(i,:))&&(norm(xtrial(i,:)-x(i,:))>0))||(any(i==id_pop_act_subpr)&&g_fun(ftrial(i,:),lambda(act_subpr(id_pop_act_subpr==i),:),z,zstar)<g_fun(f(i,:),lambda(act_subpr(id_pop_act_subpr==i),:),z,zstar)) % if agent performing local search has moved and it's objective function value has improved OR if this agent was selected to solve a sub-problem and it's partial objective function value is better than previous one
                                
                                impr(i)=1;
                                MBH_positions = [MBH_positions; xtrial(i,:)];                    %avoids cascade of fmincon usage, should be done selectively if the improvement was too little for the effort. How do I measure it?
                                
                                fprintf('Success!!!\n')
                                fprintf('Old solution = ');
                                
                                for j=1:mfit-1
                                    
                                    fprintf('%f ',f(i,j)) ;
                                    
                                end
                                
                                fprintf('%f)\n',f(i,end));
                                fprintf('New solution = ');
                                
                                for j=1:mfit-1
                                    
                                    fprintf('%f ',ftrial(i,j)) ;
                                    
                                end
                                
                                fprintf('%f)\n',ftrial(i,end));
                                
                                % delays next MBH as much as possible
                                rho(i,1)=params.rhoini;
                                rho(i,2)=0;
                                
                            end
                            
                        else
                            
                            fprintf('Fail...\n')
                            
                        end
                        
                    end
                    
                else
                    
                    % look if this position has not already been the starting point
                    % for fmincon
                    
                    if ~any(all(repmat(x(i,:),size(MBH_positions,1),1)==MBH_positions,2))
                        
                        MBH_positions = [MBH_positions; x(i,:)];                    %add this position in the list, to avoid repeating it
                        
                        fprintf('Running MBH on agent %d, lambda = (',i);
                        
                        for j =1:length(lambda(act_subpr(id_pop_act_subpr==i),:))-1
                            
                            fprintf('%f ',lambda(act_subpr(id_pop_act_subpr==i),j));
                            
                        end
                        
                        fprintf('%f)\n',lambda(act_subpr(id_pop_act_subpr==i),end));
                        niter=0;
                        fu=Inf;
                        
                        fcurrent=g_MBHfun2(x(i,:),params.func,lambda(act_subpr(id_pop_act_subpr==i),:),z,zstar,params.arg);
                        
                        while fu>=fcurrent&&niter<params.MBHflag
                            
                            niter=niter+1;
                            fprintf('Iter %d\n',niter)
                            rvlb=x(i,:)-rho(i,1)*Delta;
                            rvub=x(i,:)+rho(i,1)*Delta;
                            rvlb=max([rvlb;params.vlb]);
                            rvub=min([rvub;params.vub]);
                            
                            if niter>1
                                
                                xtrial(i,:)=x(i,:)+(2*rand-1)*rho(i,1)*Delta;
                                xtrial(i,:)=max([xtrial(i,:);rvlb]);
                                xtrial(i,:)=min([xtrial(i,:);rvub]);
                                
                            else
                                
                                xtrial(i,:) = x(i,:);
                                
                            end
                            
                            ftrial(i,:)=params.func(xtrial(i,:),params.arg{:});
                            
                            xt=xtrial(i,:);
                            
                            [u,fu,~,output]=fmincon(@(xt)g_MBHfun2(xt,params.func,lambda(act_subpr(id_pop_act_subpr==i),:),z,zstar,params.arg),xt,[],[],[],[],params.vlb,params.vub,[],foptionsNLP);
                            nfeval=nfeval+output.funcCount+1;
                            
                            if params.optimal_control==0
                                
                                xtrial(i,:) = u;
                                ftrial(i,:) = params.func(u,params.arg{:});
                                
                            else
                                
                                [ftrial(i,:),xtrial(i,:)] = params.func(u,params.arg{:});
                                
                            end
                            
                            if any(ftrial(i,:)<f(i,:))                                  % if any component of the new trial is better than the same component of the previous x, we have an improvement!
                                
                                discarded.f=[discarded.f; ftrial(i,:)];
                                discarded.x=[discarded.x; xtrial(i,:)];
                                discarded.c=[discarded.c; maxC(i)];
                                vtrial(i,:) = xtrial(i,:)-x(i,:);
                                
                                if (all(ftrial(i,:)<=f(i,:))&&(norm(xtrial(i,:)-x(i,:))>0))||(any(i==id_pop_act_subpr)&&g_fun(ftrial(i,:),lambda(act_subpr(id_pop_act_subpr==i),:),z,zstar)<g_fun(f(i,:),lambda(act_subpr(id_pop_act_subpr==i),:),z,zstar)) % if agent performing local search has moved and it's objective function value has improved OR if this agent was selected to solve a sub-problem and it's partial objective function value is better than previous one
                                    
                                    impr(i)=1;
                                    MBH_positions = [MBH_positions; xtrial(i,:)];                    %avoids cascade of fmincon usage, should be done selectively if the improvement was too little for the effort. How do I measure it?
                                    
                                    fprintf('Success!!!\n')
                                    fprintf('Old solution = ');
                                    
                                    for j=1:mfit-1
                                        
                                        fprintf('%f ',f(i,j)) ;
                                        
                                    end
                                    
                                    fprintf('%f)\n',f(i,end));
                                    fprintf('New solution = ');
                                    
                                    for j=1:mfit-1
                                        
                                        fprintf('%f ',ftrial(i,j)) ;
                                        
                                    end
                                    
                                    fprintf('%f)\n',ftrial(i,end));
                                    
                                    % delays next MBH as much as possible
                                    rho(i,1)=params.rhoini;
                                    rho(i,2)=0;
                                    
                                end
                                
                            else
                                
                                fprintf('Fail...\n')
                                
                            end
                            
                        end
                        
                    end
                    
                end
                
            end
            
        end
    end
    
    %% RADIUS CONTRACTION IF EVERYTHING FAILED, TO NARROW SEARCH AREA
    
    %     if all(params.pat_search_strategy=='standard')
    
    if impr(i)==0                                                          % if nothing worked, contract rho (seems to never kick in!)
        
        if strcmp(params.pat_search_strategy,'MADS')
            
            if rho(i,2)<params.max_rho_contr                                   % unless it has contracted 3 times, in that case restart (seems to never kick in...)
                
                rho(i,1) = rho(i,1)/4;
                rho(i,2) = rho(i,2)+1;
                
            end
            
        else
            
            rho(i,1)=rho(i,1)*params.contr_ratio;
            rho(i,2)=rho(i,2)+1;
            
            if rho(i,2)>params.max_rho_contr                                   % unless it has contracted 3 times, in that case restart (seems to never kick in...)
                
                rho(i,1)=params.rhoini;
                rho(i,2)=0;
                
            end
            
        end
        
    else
        
        %if all(ftrial(i,:)<f(i,:))
        
        if strcmp(params.pat_search_strategy,'MADS')
            
            rho(i,1)=rho(i,1)*4;
            rho(i,2)=rho(i,2)-1;
            
            if rho(i,2)<0                                               % unless it has contracted 3 times, in that case restart (seems to never kick in...)
                
                rho(i,1)=params.rhoini;
                rho(i,2)=0;
                
            end
            
            
        else
            
            rho(i,1)=rho(i,1)/params.contr_ratio;
            rho(i,2)=rho(i,2)-1;
            
            if rho(i,2)<0                                               % unless it has contracted 3 times, in that case restart (seems to never kick in...)
                
                rho(i,1)=params.rhoini;
                rho(i,2)=0;
                
            end
            
        end
        
        %end
        
    end
    
    %     end
    
    %     if all(params.pat_search_strategy=='tracking')
    %
    %         if impr(i)==0 && length(patdirs(i).avail)<1                          % if nothing worked, contract rho (seems to never kick in!)
    %
    %             rho(i,1)=rho(i,1)*params.contr_ratio;
    %             rho(i,2)=rho(i,2)+1;
    %
    %             if rho(i,2)>params.max_rho_contr                                   % unless it has contracted 3 times, in that case restart (seems to never kick in...)
    %
    %                 rho(i,1)=params.rhoini;
    %                 rho(i,2)=0;
    %
    %             end
    %
    %         else
    %
    %             if all(ftrial(i,:)<f(i,:))
    %
    %                 rho(i,1)=rho(i,1)/params.contr_ratio;
    %                 rho(i,2)=rho(i,2)-1;
    %
    %                 if rho(i,2)<0                                               % unless it has contracted 3 times, in that case restart (seems to never kick in...)
    %
    %                     rho(i,1)=params.rhoini;
    %                     rho(i,2)=0;
    %
    %                 end
    %
    %             end
    %
    %         end
    %
    %     end
    
end

return
