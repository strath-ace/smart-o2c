function [act_subpr,id_pop_act_subpr,pigr,lambda_f_best]=upd_act_subpr(pigr,lambda,lambda_f_best_old,x,f,memories,z,n_agents_subpr,maxs)
<<<<<<< HEAD:Optimisation/MACS/upd_act_subpr.m

=======
>>>>>>> 5b7361d93c9119cf1d2e9e6c885bed93f924d71b:Optimisation/MACS/upd_act_subpr.m
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
<<<<<<< HEAD:Optimisation/MACS/upd_act_subpr.m
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%
% Update subproblems (currently no longer used)

=======
%-----------Copyright (C) 2016 University of Strathclyde-------------
%
%
%
>>>>>>> 5b7361d93c9119cf1d2e9e6c885bed93f924d71b:Optimisation/MACS/upd_act_subpr.m
% This function tires to assess which subproblems are improving and selects
% new promising subproblems to investigate if some older ones are not
% progressing. It also associates the best candidate agent to solve each
% new subproblem, and updates the values of the utility function associated
% to each subproblem
%
% [act_subpr,id_pop_act_subpr,pigr,lambda_f_best,x,f]=upd_act_subpr(pigr,lambda,lambda_f_best_old,x,f,memories,z,n_agents_subpr)
%
% INPUTS
%           pigr                :       vector of values of utility function associated to each subproblem
%           lambda              :       vector of relative weights of each subproblem
%           lambda_f_best_old   :       vector of objective function values of agent best solving i-th subproblem
%           x                   :       positions of agents
%           f                   :       objective function values of agents
%           memories            :       archive
%           z                   :       min value of each objective function value found so far
%           n_agents_subpr      :       number of agents performing tackling subproblems
%
% OUTPUTS
%           act_subpr           :       updated vector of IDs of subproblems being investigated
%           id_pop_act_subpr    :       updated vector of IDs of agents tackling i-th subproblem
%           pigr                :       updated value of utility function associated to each subproblem
%           lambda_f_best       :       updated vector of objective function values of agents best solving i-th subproblem


n_lambda=size(lambda,1);                                                    % number of subproblems
lx=size(x,2);                                                               % number of pararameters
mfit=size(f,2);                                                             % number of objective functions
lambda_x_best=zeros(n_lambda,lx);                                           % position of agents which best solve each subproblem
lambda_f_best=zeros(n_lambda,mfit);                                         % objective function values of agents which best solve each subproblem

for i=1:n_lambda                                                            % for every subproblem
       
    g_tmp = hp_g_fun(memories(:,lx+1:lx+mfit),lambda(i,:),z);               % high performance evaluation of g_fun for all agents and current value of lambda
    
    [~,tmp_id]=min(g_tmp);                                                  % index of min of gfun, i.e gfun associated to agent which best solves current problem
    pigr_thres=max(g_tmp);                                                  % threshold value for pigr, equal to max gtmp (i.e gfun associated to agent with worth solution of current problem)
    lambda_x_best(i,:)=memories(tmp_id,1:lx);                               % vector of position of agents which best solve each subproblem
    lambda_f_best(i,:)=memories(tmp_id,lx+1:lx+mfit);                       % vector of objective function values of agents which best solve each subproblem
    old=g_fun(lambda_f_best_old(i,:),lambda(i,:),z);                        % compute gfun associated to old best values
    new=g_fun(lambda_f_best(i,:),lambda(i,:),z);                            % compute gfun associated to new best values
    Delta=(old-new);                                                        % variation of gfun, positive if improvement (gfun must diminish)
    
    if Delta>pigr_thres                                                     % if new gfun of best agent for this problem is lower than gfun of old best agent for this problem, and their difference is greater then the max gfun of current agents
    
        pigr(i)=1;                                                          % keep utility function to 1
        
    else                                                                    % if too little or no improvement is measured
        
        pigr(i)=(0.95+0.05*Delta/pigr_thres)*pigr(i);                       % utility function diminishes
        
    end
    
end

% Selects the new population of subproblems to be treated
if mfit>1                                                                   % if problem is multiobjective

    tourn_size=round(n_lambda/60);                                          % why??????
    
else                                                                        % if problem is single objective
    
    tourn_size=round(n_lambda/2);                                           % why??????
    
end

tmp_id=(mfit+1):n_lambda;                                                   % vector of subproblems IDs, excluding first orthogonal ones
act_subpr=[1:mfit zeros(1,n_agents_subpr-mfit)];                            % vector of next active subproblems, for sure composed by mfit orthogonal subproblems

for i=mfit+1:n_agents_subpr                                                 % for all subproblems (remember, only one agent per subproblem)

    sel=randperm(length(tmp_id));                                           % shuffle the vector of (remaining) non orthogonal subproblems
    sel=sel(1:min([tourn_size length(sel)]));                               % pick no more than tourn_size of non orthogonal subproblems (or pick them all if they are less than tourn size)
    [~,id]=max(pigr(tmp_id(sel)));                                          % within the above mentioned selection, pick the ID of the one with highest utility function 
    %id = 1;
    act_subpr(i)=tmp_id(sel(id));                                           % the i-th active subproblem is the one just chosen (i.e the one with the highest utility function)
    tmp_id=tmp_id(tmp_id~=act_subpr(i));                                    % remove the subproblem just selected from the list of possible subproblems
    
end

id_pop_act_subpr=zeros(1,n_agents_subpr);                                   % initialize the vector containing the IDs of agents which will explore subproblems (one each)
list=1:size(f,1);                                                           % create a list of IDs from 1 to num of agents

for i=1:n_agents_subpr                                                      % for every agent tackling a subproblem

    g_tmp=Inf;                                                              % g is initialized as Inf
      
    for j=list                                                              % for every fun value of current agents
    
        tmp=g_fun(f(j,:),lambda(act_subpr(i),:),z,maxs);                         % evaluate it's g fun associated to this subproblem
        
        if tmp<g_tmp                                                        % if it's better than the old one
        
            g_tmp=tmp;                                                      % substitute it
            id_pop_act_subpr(i)=j;                                          % and associate current j-th element to i-th subproblem
            
        end
        
    end
    
    list=list(list~=id_pop_act_subpr(i));                                   % remove current j-th agent from the list of agents that is possible to associate to subproblems
    
end

return
