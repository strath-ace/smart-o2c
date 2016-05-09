function [Solutions, ListNodes, Agents, agentdeathflag] = AgentMovement(Inputs, Solutions, ListNodes, Agents, agent)
% This function is used to move agents through the graph/tree
%
% Inputs:
% * Inputs      : Structure containing the PhysarumSolver inputs
% * ListNodes       : Structure containing the graph
% * Agents      : Structure containing the Agents
% * agent       : Cell with the current agents' name
% 
% Outputs:
% * ListNodes       : Structure containing the graph
% * Agents      : The agents with their new positions
% * agentdeathflag : Flag that shows whether the agent has 'died'
%
% Author: Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk

%Initialize the agent death flag
agentdeathflag = 0;

%For ease of reading, define the current agent and the current node as
%variables
currentagent = char(agent);
currentnode = char(Agents.(currentagent).currentNode);
disp(strcat([datestr(now),' === Moved to',' ',currentnode]));

if ((isempty(ListNodes.(currentnode).possibledecisions)) || (sum(ListNodes.(currentnode).VisitsLeft) == 0))
    
    %If there are no more possible decisions or visists left, set the death
    %flag to 1
    agentdeathflag = 1;
    
    %Save the solution
    Solutions.Nodes = [Solutions.Nodes; {[Agents.(agent).previousListNodes {Agents.(agent).currentNode}]}]';
    return
end

%Generate a random number
p = rand;

%check whether the current node has children & whether the random number
%falls outside of the probability margin
if (p>Inputs.RamificationProbability && ~isempty(ListNodes.(currentnode).children))
    
    %Pre-allocate and calculate the probabilities to transverse to each node. The smaller
    %the cost (higher flux), the higher the probability
    childfluxes = zeros(1,length(ListNodes.(currentnode).children));
    for i = 1:length(ListNodes.(currentnode).children)
        childfluxes(i) = ListNodes.(char(ListNodes.(currentnode).children(i))).flux;
    end
    problist = childfluxes./(sum(childfluxes));
    problist = problist./sum(problist);

    %Choose one of the node's children. Cell structure is used to
    %circumvent issues with node selection if only 1 child is present and
    %hocsenindex = 1 (it will otherwise only return the first letter)
    chosenindex = find(rand<cumsum(problist),1,'first');
    nodechildren = cell(ListNodes.(currentnode).children);
    chosennode = nodechildren{chosenindex};
    
    %Move agent to the new node
    Agents.(currentagent).previousListNodes = [Agents.(currentagent).previousListNodes {currentnode}];
    Agents.(currentagent).currentNode = chosennode;
    Agents.(currentagent).previouscosts = [Agents.(currentagent).previouscosts ListNodes.(currentnode).length];
else
    
    %If there are no children or if the random number falls within the
    %probability margin, ramificate towards a new node
    [ListNodes, Solutions, Agents, agentdeathflag] = Ramification(Inputs, Solutions, ListNodes, Agents, currentagent);
end
   