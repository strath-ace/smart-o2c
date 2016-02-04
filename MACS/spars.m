function A = spars (A,x,f,c,z)

% Substitutes agents x in archive A if doing so would reduce clustering of 
% agents in archive
%
% A = spars (A,x,f,c,z)
%
% INPUTS
%           A       :       archive
%           x       :       position of agents
%           f       :       objective function values of agents
%           c       :       constraint violations of agents
%           z       :       minimas of objective functions
%
% OUTPUTS
%           A       :       updated archive


[nX,lx] = size(x);                                                          % number of agents, and number of paremeters for each agent

for i = 1:nX                                                                % for each agent
    
    if all(f(i,:)>z)                                                        % if all objective values of i-th object are greater than current minima (i.e it's dominant but NOT optimal wrt some subproblem). Local extrema should have already been added so skip them

        [DeltaX,sD]=deltacomp(A,x(i,:),0);                                  % compute and sort Euclidean distances from i-th agent to all agents in memory
        [DeltaA]=deltacomp(A,A(sD(1),1:lx),sD(1));                          % compute and sort Euclidean distances from element closest to i-th agent (excluding i-th agent itself) to all agents in memory
        
        if DeltaX(2)>DeltaA(1)                                              % if distance from this agent to second nearest agent in memory is greater than distance from closest agent in memory to second closest agent in memory
        
            A(sD(1),:)=[x(i,:) f(i,:) 0 c(i)];                              % substitute this agents to agent in memory closest to it (i.e, this move reduces clustering)
            
        end
        
    end
    
end

return