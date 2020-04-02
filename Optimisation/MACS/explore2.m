function [xtrial,vtrial,ftrial,maxC,nfeval,discarded,rho,patdirs,MBH_positions,MADS_dirs]=explore2(memories,x,v,f,cid,nfeval,lambda,act_subpr,id_pop_act_subpr,z,zstar,rho,patdirs,pigr,MBH_positions,MADS_dirs,local_only,int_loc_opt,params)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%
% Performs individual moves

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
foptionsNLP=params.mbh_options;

%% MAIN LOOP
impr = zeros(n_a,1);

%fprintf(['Performing individual moves: \n']);

for i=1:n_a                                                                 % for all agents which will perform LOCAL ACTIONS
    
    
    fprintf(['Performing individual moves on agent ' ,num2str(i),',']);
    %fprintf(['(', num2str(rho(i,2)) , '):\n ']);
    
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
            
            fprintf('Performing intertia: ');
            
            %% EVALUATE MOVE
            
            if params.cp==0                                                 % if unconstrained problem
                
                if params.bilevel==1
                    
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
                        fprintf('V\n');     % found dominating point or point with better scalarisation
                        
                    else
                        
                        fprintf('v\n');     % found non dominated point (wrt this agent only)
                        
                    end
                    
                    % IDEA, BETTER TO CONTINUE THE LOOP THAN CHECK THIS VAR ONWARDS, SAVES TIME AND READABILITY
                    %else
                    
                    %vtrial(i,:) = 0*vtrial(i,:);
                    
                else
                    
                    fprintf('X\n');         % found dominated point
                    
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
        
        fprintf(['Performing Pattern Search (up to ',num2str(n_ids),' coordinates): ']);
        
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
                    
                    if params.bilevel ==1
                        
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
                            fprintf('V\n');
                            
                        else
                            
                            fprintf('v');
                            
                        end
                        
                    else
                        
                        fprintf('X');
                        
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
                        
                        if params.bilevel ==1
                            
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
                                fprintf('V\n');
                                
                            else
                                
                                fprintf('v,');
                                
                            end
                            
                        else
                            
                            fprintf('X,');
                            
                        end
                        
                        if j==n_ids
                            
                            fprintf('\n');
                            
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
                    
                    if params.bilevel ==1
                        
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
                    
                    if params.bilevel ==1
                        
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
                    
                    if params.bilevel ==1
                        
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
        
        fprintf('Performing Differential Evolution: ');
        
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
                
                if params.bilevel ==1
                    
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
                        fprintf('V\n');
                        
                        if params.v==1
                            
                            vtrial(i,:) = xtrial(i,:)-x(i,:);
                            
                        end
                        
                    else
                        
                        fprintf('v\n');
                        
                    end
                    
                else
                    
                    fprintf('X\n');
                    
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
    fprintf('\n');
    
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
    
    %% MBH
    
    if impr(i)==0
        
        if params.optimal_control==1
            
            if any(i==id_pop_act_subpr)&& params.MBHflag>0 && local_only%any(i==id_pop_act_subpr)&& params.MBHflag>0 && (rho(i,2)==params.max_rho_contr || local_only) %&&pigr(act_subpr(id_pop_act_subpr==i))<0.8%<- TO BE REVISED WITH A PROPPER PARAMETER TO BE SET BY THE USER
                
                % look if this position has not already been the starting point
                % for fmincon
                
                fprintf('Running MBH on agent %d, lambda = (',i);
                
                for j =1:length(lambda(act_subpr(id_pop_act_subpr==i),:))-1
                
                    fprintf('%f ',lambda(act_subpr(id_pop_act_subpr==i),j));
                
                end
                
                fprintf('%f)\n',lambda(act_subpr(id_pop_act_subpr==i),end));
                niter=0;
                fu=Inf;
                rand_restart = 0;
                change_dir = 0;
                % look if current starting position has already been used
                
                selvec = all(repmat(x(i,:),size(MBH_positions,1),1)==MBH_positions(:,1:lx),2);
                
                if any(selvec)
                    
                    % if this position was used by an orht-following agent
                    if act_subpr(id_pop_act_subpr==i)<=mfit
                        
                        thislist = MBH_positions(selvec,:);
                        
                        % check which direction was this position used for
                        % and if it matches with this agent's direction
                        
                        thislist_dir = thislist(:,lx+1:lx+mfit);
                        availdirs = [lambda(act_subpr(id_pop_act_subpr==i),:); 1-lambda(act_subpr(id_pop_act_subpr==i),:)];
                        availdirs(2,act_subpr(id_pop_act_subpr==i)) = -1;
                        
                        if size(intersect(thislist_dir,availdirs,'rows','stable'),1)==1
                            
                            change_dir = 1;
                            %fprintf('This agent already failed the local search from its starting position... trying another direction\n')
                            
                        else
                            
                            rand_restart = 1;
                            %fprintf('This agent already failed the local search from its starting position... restarting in the local neighbourhood\n')
                            
                        end
                        
                    else
                        
                        %non orth agents have no possibility to change
                        %direction (yet), so just regenerate them
                        rand_restart = 1;
                        %fprintf('This agent already failed the local search from its starting position... restarting in the local neighbourhood\n')
                        
                    end
                    
                end
                
                fcurrent=1;
                
                if ~rand_restart && mod(int_loc_opt,2)==0
                    
                    change_dir = 1;
                    
                end
                
                while fu>=fcurrent&&niter<params.MBHflag
                    
                    niter=niter+1;
                    %fprintf('Iter %d\n',niter)
                    rvlb=x(i,:);
                    rvub=x(i,:);
                    rvlb(params.vars_to_opt) = x(i,params.vars_to_opt)-rho(i,1)*Delta(params.vars_to_opt);
                    rvub(params.vars_to_opt) = x(i,params.vars_to_opt)+rho(i,1)*Delta(params.vars_to_opt);
                    
                    rvlb=max([rvlb;params.vlb]);
                    rvub=min([rvub;params.vub]);
                    
                    if (niter>1) || (rand_restart==1)
                        
                        xtrial(i,:)=x(i,:);
                        
                        xtrial(i,params.vars_to_opt) = xtrial(i,params.vars_to_opt)+(2*rand-1)*rho(i,1)*Delta(params.vars_to_opt);
                        
                        xtrial(i,:)=max([xtrial(i,:);rvlb]);
                        xtrial(i,:)=min([xtrial(i,:);rvub]);
                        
                        if params.bilevel ==1
                            
                            [~,xtrial(i,:)]=params.func(xtrial(i,:),params.arg{:});         % eval f on given position
                            nfeval = nfeval+1;
                            
                        end
                        
                    else
                        
                        xtrial(i,:) = x(i,:);
                        
                    end
                    
                    % avoid starting on the bounds for interior-point
                    
                    if strcmp(foptionsNLP.Algorithm,'interior-point')
                        
                        xtrial(i,xtrial(i,:)==params.vlb)= params.vlb(xtrial(i,:)==params.vlb)+1e-6;
                        xtrial(i,xtrial(i,:)==params.vub)= params.vub(xtrial(i,:)==params.vub)-1e-6;
                        
                    end
                    
                    if act_subpr(id_pop_act_subpr==i)<=mfit
                        
                        if change_dir==0
                            
                            % pursue orht subproblem
                            zz = f(i,:)-2*(zstar-z);
                            ll = lambda(act_subpr(id_pop_act_subpr==i),:);
                            
                        else
                            
                            % try to minimise other objectives while
                            % keeping the value of its primary obj
                            
                            zz = z;
                            ll = 1-lambda(act_subpr(id_pop_act_subpr==i),:);
                            ll(act_subpr(id_pop_act_subpr==i)) = -1;
                            
                        end
                        
                    else
                        
                        zz = f(i,:)-2*(zstar-z);
                        ll = ones(1,mfit)./norm(ones(1,mfit));
                        
                    end
                    
                    % xtrial needs normalisation
                    
                    xstart = xtrial(i,:);
                    normx = xtrial(i,:)./params.oc.problem.scales.scale_opt';%normalise_vect(xtrial(i,:),params.oc.scale_opt,params.oc.offs_opt);
                    xt=[max(ll) normx];
                    
                    fcurrent=max(ll);
                    zzstar = f(i,:);    % improvement is always wrt ORIGINAL position, not restarted position                    
                    
                    try
                        % for some mysterious reason, i need a bogus command before calling fmincon or the try statement won't even trigger...
                        
                        tstart = tic;                    
                        timeout = 0;
                        [u,fu,~,output]=fmincon(@smooth_cheb_fun,xt,[],[],[],[],[0 params.oc.problem.norm_lb'],[max(ll)+0.1 params.oc.problem.norm_ub'],@(xt) params.smooth_scal_constr_fun(xt,ll,zz,zzstar,params,tstart),foptionsNLP);
                        
                    catch EX
                        
                        if strcmp(EX.identifier,'multiphase_constr_full:timeOutReached')
                            
                            timeout = 1;
                            
                        else
                           
                            rethrow(EX);
                            
                        end
                        
                    end
                    
                    if ~timeout
                        
                        output.message
                        nfeval=nfeval+output.funcCount+1;
                        denormx = u(2:end).*params.oc.problem.scales.scale_opt';%denormalise_vect(u(2:end),params.oc.scale_opt,params.oc.offs_opt);
                        
                        if params.bilevel==0
                            
                            xtrial(i,:) = denormx;
                            ftrial(i,:) = params.func(denormx,params.arg{:});
                            
                        else
                            
                            [ftrial(i,:),xtrial(i,:)] = params.func(denormx,params.arg{:}); % eval objectives of this solution (not possible to retrieve them from fu)
                            nfeval = nfeval+1;
                            
                        end
                        
                        %fprintf('Old solution = ');
                        
                        %for j=1:mfit-1
                        
                        %fprintf('%f ',f(i,j)) ;
                        
                        %end
                        
                        %fprintf('%f)\n',f(i,end));
                        %fprintf('New solution = ');
                        
                        %for j=1:mfit-1
                        
                        %    fprintf('%f ',ftrial(i,j)) ;
                        
                        %end
                        
                        %fprintf('%f)\n',ftrial(i,end));
                        
                        if any(ftrial(i,:)<f(i,:))                                  % if any component of the new trial is better than the same component of the previous x, we have an improvement!
                            
                            discarded.f=[discarded.f; ftrial(i,:)];
                            discarded.x=[discarded.x; xtrial(i,:)];
                            discarded.c=[discarded.c; maxC(i)];
                            
                            if (all(ftrial(i,:)<=f(i,:))&&(norm(xtrial(i,:)-x(i,:))>0))||(any(i==id_pop_act_subpr)&&g_fun(ftrial(i,:),lambda(act_subpr(id_pop_act_subpr==i),:),z,zstar)<g_fun(f(i,:),lambda(act_subpr(id_pop_act_subpr==i),:),z,zstar)) % if agent performing local search has moved and it's objective function value has improved OR if this agent was selected to solve a sub-problem and it's partial objective function value is better than previous one
                                
                                impr(i)=1;
                                vtrial(i,:) = xtrial(i,:)-x(i,:);
                                %fprintf('Success!!!\n')
                                
                            else
                                
                                %fprintf('MBH found a non dominated solution in a direction different than the one associated to the agent...\n');
                                
                            end
                            
                        else
                            
                            %fprintf('Fail...\n')
                            
                        end
                        
                        
                    end
                    
                    MBH_positions = [MBH_positions; xstart ll];
                    
                end
                
            end
            
        else
            
            % check all code below!!!
            
            if params.MBHflag>0 && local_only%any(i==id_pop_act_subpr)&& params.MBHflag>0 && (rho(i,2)==params.max_rho_contr || local_only) %&&pigr(act_subpr(id_pop_act_subpr==i))<0.8%<- TO BE REVISED WITH A PROPPER PARAMETER TO BE SET BY THE USER
                
                % look if this position has not already been the starting point
                
                %fprintf('Running MBH on agent %d\n)',i);
                
                niter=0;
                fu=Inf;
                rand_restart = 0;
                change_dir = 0;
                % look if current starting position has already been used
                
                selvec = all(repmat(x(i,:),size(MBH_positions,1),1)==MBH_positions(:,1:lx),2);
                
                if any(selvec) || any(f(i,:)==z)
                    
                    % if this is a social agent
                    
                    if any(id_pop_act_subpr==i)
                        
                        % if this position was used by an orht-following agent
                        
                        if act_subpr(id_pop_act_subpr==i)<=mfit
                            
                            thislist = MBH_positions(selvec,:);
                            
                            % check which direction was this position used for
                            % and if it matches with this agent's direction
                            
                            thislist_dir = thislist(:,lx+1:lx+mfit);
                            availdirs = [lambda(act_subpr(id_pop_act_subpr==i),:); 1-lambda(act_subpr(id_pop_act_subpr==i),:)];
                            availdirs(2,act_subpr(id_pop_act_subpr==i)) = -1;
                            
                            if size(intersect(thislist_dir,availdirs,'rows','stable'),1)==1
                                
                                change_dir = 1;
                                %fprintf('This agent already failed the local search from its starting position... trying another direction\n')
                                
                            else
                                
                                rand_restart = 1;
                                %fprintf('This agent already failed the local search from its starting position... restarting in the local neighbourhood\n')
                                
                            end
                            
                        else
                            
                            % non orth agents have no possibility to change
                            % direction (yet), so just regenerate them
                            rand_restart = 1;
                            %fprintf('This agent already failed the local search from its starting position... restarting in the local neighbourhood\n')
                            
                        end
                        
                    else
                        
                        % non social agents have no possibility to change
                        % direction , so just regenerate them
                        rand_restart = 1;
                        %fprintf('This agent already failed the local search from its starting position... restarting in the local neighbourhood\n')
                        
                    end
                    
                end
                
                if ~rand_restart && mod(int_loc_opt,2)==0
                    
                    change_dir = 1;
                    
                end
                
                fcurrent=1;
                
                while fu>=fcurrent&&niter<params.MBHflag
                    
                    niter=niter+1;
                    %fprintf('Iter %d\n',niter)
                    rvlb=x(i,:);
                    rvub=x(i,:);
                    rvlb(params.vars_to_opt) = x(i,params.vars_to_opt)-rho(i,1)*Delta(params.vars_to_opt);
                    rvub(params.vars_to_opt) = x(i,params.vars_to_opt)+rho(i,1)*Delta(params.vars_to_opt);
                    
                    rvlb=max([rvlb;params.vlb]);
                    rvub=min([rvub;params.vub]);
                    
                    if (niter>1) || (rand_restart==1)
                        
                        xtrial(i,:)=x(i,:);
                        
                        xtrial(i,params.vars_to_opt) = xtrial(i,params.vars_to_opt)+(2*rand(size(Delta(params.vars_to_opt)))-1).*rho(i,1).*Delta(params.vars_to_opt);
                        
                        xtrial(i,:)=max([xtrial(i,:);rvlb]);
                        xtrial(i,:)=min([xtrial(i,:);rvub]);
                        
                        if params.bilevel ==1
                            
                            [~,xtrial(i,:)]=params.func(xtrial(i,:),params.arg{:});         % eval f on given position
                            nfeval = nfeval+1;
                            
                        end
                        
                        
                    else
                        
                        xtrial(i,:) = x(i,:);
                        
                    end
                    
                    if any(id_pop_act_subpr==i)
                        
                        if act_subpr(id_pop_act_subpr==i)<=mfit
                            
                            if change_dir==0
                                
                                % pursue orht subproblem
                                zz = f(i,:)-2*(zstar-z);
                                ll = lambda(act_subpr(id_pop_act_subpr==i),:);
                                
                            else
                                
                                % try to minimise other objectives while
                                % keeping the value of its primary obj
                                
                                zz = z;
                                ll = 1-lambda(act_subpr(id_pop_act_subpr==i),:);
                                ll(act_subpr(id_pop_act_subpr==i)) = -1;
                                
                            end
                            
                        else
                            
                            zz = f(i,:)-2*(zstar-z);
                            ll = ones(1,mfit)./norm(ones(1,mfit));
                            
                        end
                        
                    else
                        
                        zz = f(i,:)-2*(zstar-z);
                        ll = ones(1,mfit)./norm(ones(1,mfit));
                        
                    end
                    
                    % xtrial needs normalisation
                    xstart = xtrial(i,:);
                    normx = xtrial(i,:)./params.scales;
                    
                    xt=[max(ll) normx];
                    
                    fcurrent=max(ll);
                    zzstar = f(i,:);
                    
                    timeout = 0;
                    
                    [u,fu,~,output]=fmincon(@smooth_cheb_fun,xt,[],[],[],[],[0 params.norm_vlb],[max(ll)+0.1 params.norm_vub],@(xt) params.smooth_scal_constr_fun(xt,params.mbh_func,params.mbh_cfunc,ll,zz,zzstar,params.scales,params.arg{:}),foptionsNLP);
                    
                    output.message
                    nfeval=nfeval+output.funcCount+1;
                    
                    xtrial(i,:) = u(2:end).*params.scales;%denormalise_vect(u(2:end),params.oc.scale_opt,params.oc.offs_opt);
                    
                    % hard clipping (might get out of bounds by eps, and
                    % still break some functions)
                    
                    xtrial(i,:)=max([xtrial(i,:);params.vlb]);
                    xtrial(i,:)=min([xtrial(i,:);params.vub]);
                    
                    
                    if params.bilevel==0
                        
                        ftrial(i,:) = params.func(xtrial(i,:),params.arg{:});
                        
                    else
                        
                        [ftrial(i,:),xtrial(i,:)] = params.func(xtrial(i,:),params.arg{:}); % eval objectives of this solution (not possible to retrieve them from fu)
                        nfeval = nfeval+1;
                        
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
                        
                        if (all(ftrial(i,:)<=f(i,:)))||(any(i==id_pop_act_subpr)&&g_fun(ftrial(i,:),lambda(act_subpr(id_pop_act_subpr==i),:),z,zstar)<g_fun(f(i,:),lambda(act_subpr(id_pop_act_subpr==i),:),z,zstar)) % if agent performing local search has moved and it's objective function value has improved OR if this agent was selected to solve a sub-problem and it's partial objective function value is better than previous one
                            
                            impr(i)=1;
                            vtrial(i,:) = xtrial(i,:)-x(i,:);
                            %fprintf('Success!!!\n')
                            
                        else
                            
                            %fprintf('MBH found a non dominated solution in a direction different than the one associated to the agent...\n');
                            
                        end
                        
                    else
                        
                        %fprintf('Fail...\n')
                        
                    end
                    
                    MBH_positions = [MBH_positions; xstart ll];
                    
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
    
    %     figure()
    %     plot(f(:,1),f(:,2),'go')
    %     hold on
    %     plot(f(i,1),f(i,2),'r+')
    %     plot(discarded.f(ndisc+1:end,1),discarded.f(ndisc+1:end,2),'b.')
    %     hold off
    %     ndisc = size(discarded.f,1);
    %keyboard
    
    
end


%fprintf('\n');

return
