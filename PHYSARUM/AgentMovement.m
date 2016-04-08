function [Nodes,Agents,agentdeathflag] = AgentMovement(Inputs,Nodes,Agents,agent)
% This function is used to move agents through the graph/tree
%
% Inputs:
% * Inputs      : Structure containing the PhysarumSolver inputs
% * Nodes       : Structure containing the graph
% * Agents      : Structure containing the Agents & their characteristics
% * agent       : Cell with the current agents' name
% 
% Outputs:
% * Nodes       : Structure containing the graph
% * Agents      : The agents with their new positions
%
% Author: Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk

%Initialize the agent death flag
 agentdeathflag = 0;

%For ease of reading, define the current agent and the current node as
%variables
currentagent = char(agent);
currentnode = char(Agents.(currentagent).currentNode);

%Generate a random number
p = rand;

%check whether the current node has children & whether the random number
%falls outside of the probability margin
if (p>Inputs.RamificationProbability && ~isempty(Nodes.(currentnode).children))
    
    %Calculate the probabilities to transverse to each node. The smaller
    %the cost, the higher the probability
    problist = 1./(Nodes.(currentnode).lengths);
    problist = problist./sum(problist);

    %Choose one of the node's children. Cell structure is used to
    %circumvent issues with node selection if only 1 child is present and
    %hocsenindex = 1 (it will otherwise only return the first letter)
    chosenindex = find(rand<cumsum(problist),1,'first');
    nodechildren = cell(Nodes.(currentnode).children);
    chosennode = nodechildren{chosenindex};
    
    %Move agent to the new node
    Agents.(currentagent).previousNodes = [Agents.(currentagent).previousNodes {currentnode}];
    Agents.(currentagent).currentNode = chosennode;
    Agents.(currentagent).previouscosts = [Agents.(currentagent).previouscosts Nodes.(currentnode).lengths(chosenindex)];
else
    
    %If there are no children or if the random number falls within the
    %probability margin, ramificate towards a new node
    [Nodes,Agents,agentdeathflag] = Ramification(Inputs,Nodes,Agents,currentagent);
end
   