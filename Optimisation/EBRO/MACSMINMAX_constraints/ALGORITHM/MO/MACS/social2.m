function [xsamp,fsamp,maxC,nfeval]=social2(x,f,memories,x_DE,cid,nfeval,energy,ener2,params)

% [xsamp,fsamp,maxC,nfeval]=social(x,f,memories,x_DE,cid,nfeval,params)
%
%  INPUT
%           x                   : positions of the agents on which DE will be performed
%           f                   : objective values associated to agents on which DE will be performed
%           memories            : archive of previous solutions 
%           x_DE                : position of all agents
%           cid                 : max constraint violation for agents on which DE will be performed
%           nfeval              : number of function evaluations
%
%           params  : structure of parameters composed by
%               'F'             : deviation factor for DE
%               'CR'            : crossover factor for DE
%               'vlb'           : lower boundary
%               'vub'           : upper boundary
%               'func'          : function handle to objective functions
%               'arg'           : argument of func
%               'cp'            : flag for unconstrained/constrained problem
%               'T'             : min desired number of agents composing local network
%               'DE_strategy'   : strategy to use
%
%  OUTPUT
%           xsamp               : sampled values
%           fsamp               : objective functions associated to xsamp
%           maxC                : max constraint violation associated to xsamp
%           nfeval              : updated number of function evaluations
%
% CHANGE LOG
% Created by: Massimiliano Vasile 2010
% Revised, cleaned and optimized by: Lorenzo A. Ricciardi 2015

%% DE INITIALIZATION

params.T=max([params.T 4]);                                                 % local network size, should be at least 4 agents 
xsamp=x;                                                                    % position of agents on which DE will act
fsamp=f;                                                                    % objective function value of agents on which DE will act
csamp=zeros(size(x,1),1);                                                   % constraint violation of agents on which DE will act
maxC=csamp;                                                                 % max constraint violation of agents on which DE will act
lx = length(params.vlb);                                                    % number of parameters

%memories(all(memories(:,1:lx)==repmat(x,size(memories,1),1),2),:) = [];     % consider only agents in memory whose position is different than xsamp (cool one-liner, huh? :P)
%memories = unique(memories,'rows');                                         % consider only unique agents from memory

x_DE(all(x_DE==repmat(x,size(x_DE,1),1),2),:) = [];                         % consider only agents in current population whose position is different than xsamp (cool one-liner, huh? :P)
x_DE = unique(x_DE,'rows');                                                 % consider only unique agents in current population

n_mem=size(memories,1);                                                     % number of elements in memory
pop_DE=size(x_DE,1);                                                        % number of elements from  all DE population

p_arch_vs_pop=1-exp(-n_mem/pop_DE);                                         % probability of choosing an element from the archive rather than from current DE population
arch_vs_pop=rand<p_arch_vs_pop;                                             % either 0 or 1, random. Is 1 if rand<p_arch_vs_pop


if arch_vs_pop&&n_mem>3                                                    % if we are taking samples from the archive and the archive contains at least 3 elements

    tmp_dist=(sum(((memories(:,params.id_vars_to_opt)-repmat(xsamp(params.id_vars_to_opt),n_mem,1))./repmat(params.vub(params.id_vars_to_opt)-params.vlb(params.id_vars_to_opt),n_mem,1)).^2,2)).^0.5; % fast one-liner to compute distances from all agents in memory to xsample
    [~,P]=sort(tmp_dist);                                                   % sort elements in tmp_dist, storing their position in P (THIS IS THE LOCAL SEARCH NETWORK)
    %[~,P] = sort(ener2);
    P=P(randperm(min([params.T n_mem])))';                                  % pick at least 4 elements (because T is at least 4; if not possible pick as many as tmp allows) at random from THE LOCAL SEARCH NETWORK (i.e consider the closest 4/length(tmp_dist) elements from xsamp)
    
else                                                                        % if we are not taking samples from the archive or we don't have at least 3 samples from the archive
    
    tmp_dist=(sum(((x_DE(:,params.id_vars_to_opt)-repmat(xsamp(params.id_vars_to_opt),pop_DE,1))./repmat(params.vub(params.id_vars_to_opt)-params.vlb(params.id_vars_to_opt),pop_DE,1)).^2,2)).^0.5; % fast one-liner to compute distances from all agents in population to xsample
    [~,id_pop]=sort(tmp_dist);                                              % sort tmp_dist (THIS IS THE LOCAL SEARCH NETWORK)
    id_pop=id_pop(randperm(min([params.T pop_DE])));                        % pick at least 4 elements (because T is at least 4; if not possible pick as many as tmp_dist allows) at random from THE LOCAL SEARCH NETWORK (i.e consider the closest 4/length(tmp_dist) elements from xsamp)
    
end

%% DE

switch params.DE_strategy
    
    case 'DE/rand/1'
        
        e=rand(1,lx)<params.CR;                                             % mutation vector, with probability CR for each element
    
        if (arch_vs_pop) && (n_mem>2)                                       % if we are taking samples from archives, and there are at least 3 samples in it (remember that we have excluded x itself!)
        
            xsamp(params.id_vars_to_opt)=xsamp(params.id_vars_to_opt).*(1-e)+e.*(memories(P(3),params.id_vars_to_opt)+params.F*(memories(P(1),params.id_vars_to_opt)-memories(P(2),params.id_vars_to_opt))); % xsamp is evolved with differentiation vector given by agents P(1), P(2) and P(3) FROM THE ARCHIVE
            
        else                                                                % if we are not taking samples from the archive, or it doesn't contain enough samples
            
            if size(id_pop,1)<2
                
                xsamp = xsamp;
                
            else
                if size(id_pop,1)==2
                    
                    xsamp(params.id_vars_to_opt)=xsamp(params.id_vars_to_opt).*(1-e)+e.*(params.F*(x_DE(id_pop(1),params.id_vars_to_opt)-x_DE(id_pop(2),params.id_vars_to_opt)));  % xsamp is evolved with differentiation vector given by agents id_pop(1), id_pop(2) and id_pop(3) FROM THE POPULATION
                    
                else
                    
                    xsamp(params.id_vars_to_opt)=xsamp(params.id_vars_to_opt).*(1-e)+e.*(x_DE(id_pop(3),params.id_vars_to_opt)+params.F*(x_DE(id_pop(1),params.id_vars_to_opt)-x_DE(id_pop(2),params.id_vars_to_opt)));  % xsamp is evolved with differentiation vector given by agents id_pop(1), id_pop(2) and id_pop(3) FROM THE POPULATION
                    
                end
                
            end
            
        end
        
    case 'DE/current-to-rand/1'
        
        K= rand;                                                             % create a random number K
        K2 = K;%rand;
        if (arch_vs_pop) && (n_mem>3)                                       % if we are sampling from the archive and it contains at least 3 elements
        
            xsamp(params.id_vars_to_opt)=xsamp(params.id_vars_to_opt)+K*(memories(P(3),params.id_vars_to_opt)-xsamp((params.id_vars_to_opt)))+K2*params.F*(memories(P(1),params.id_vars_to_opt)-memories(P(2),params.id_vars_to_opt)); % xsamp is evolved with differentiation vector given by agents P(1) P(2) and P(3) FROM THE ARCHIVE
            
        else
            
            if size(id_pop,1)<2
                
                xsamp = xsamp;
                
            else
                
                if size(id_pop,1)==2
                    
                    xsamp(params.id_vars_to_opt)=xsamp(params.id_vars_to_opt)+K*params.F*(x_DE(id_pop(1),params.id_vars_to_opt)-x_DE(id_pop(2),params.id_vars_to_opt));
                    
                else
                    
                    xsamp(params.id_vars_to_opt)=xsamp(params.id_vars_to_opt)+K*(x_DE(id_pop(3),params.id_vars_to_opt)-xsamp(params.id_vars_to_opt))+K2*params.F*(x_DE(id_pop(1),params.id_vars_to_opt)-x_DE(id_pop(2),params.id_vars_to_opt));
                    
                end
                
            end
            
        end
        
end

%% BOUNDS CHECK WITH RANDOM RESAMPLING IF OUT OF BOUNDS

mask = (xsamp(params.id_vars_to_opt)<params.vlb(params.id_vars_to_opt))+(xsamp(params.id_vars_to_opt)>params.vub(params.id_vars_to_opt));                                             % mask vector, contains 1 where xsamp is out of bounds
xnew = mask.*(rand(1,length(params.id_vars_to_opt)).*(params.vub(params.id_vars_to_opt)-params.vlb(params.id_vars_to_opt))+params.vlb(params.id_vars_to_opt));                                     % xnew, vector with random components where xsamp is out of bounds AND ZERO ELSEWHERE
xsamp(params.id_vars_to_opt) = xsamp(params.id_vars_to_opt).*~mask+xnew.*mask;                                            % in bounds xsamp

%% EVALUATION OF OBJECTIVE FUNCTION

if norm(xsamp-x)>0                                                          % if xsamp is different than x
    
    if params.cp==0                                                                % if problem is unconstrained
        
        if params.optimal_control ==1
        
            [fsamp,xcorr]=params.func(xsamp,params.arg{:});                                           % evaluate f
            xsamp = xcorr;
            maxC=0;
        
        else
            
            fsamp=params.func(xsamp,params.arg{:});                                           % evaluate f
            maxC=0;
            
        end
        
    else                                                                    % if it's constrained
        
        [fsamp,csamp]=params.func(xsamp,params.arg{:});                                   % evaluate f
        maxC=max(csamp);                                                    % and max constraint violation
        
        if cid<=0&&maxC>0                                                   % if previous solution was in feasible region and current one is not
            
            fsamp=f+maxC;                                                   % penalize this solution adding max constraint violation
            
        end
        
    end
    
    nfeval=nfeval+1;
    
end

%%

return