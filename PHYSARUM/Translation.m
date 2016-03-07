function [Nodes] = Translation(Inputs,Nodes,Agents)
% This function is used to move agents through the graph/tree
%
% Inputs:
% * Nodes       : Structure containing the graph
% * Inputs      : Structure containing the PhysarumSolver inputs
% * Agents      : Structure containing the agents' characteristics
%
% Outputs:
% * Nodes       : Structure containing the graph

 
ListOfAgents = fields(Agents);
for i = 1:size(ListOfAgents,1)
    
    %For easy of reading, define the name of the agent's current node
    %as a separate variable
    AgentNodeName =  char(Agents.(char(ListOfAgents(i))).currentNode);
    
    %Add agent's current node to its history of previous nodes
    Agents.(char(ListOfAgents(i))).previousNodes = [Agents.(char(ListOfAgents(i))).previousNodes {Agents.(char(ListOfAgents(i))).currentNode}];
    
    %Calculate the probabilities for the children using their fluxes
    Nodes.ListNodes.(AgentNodeName).probabilities = Nodes.ListNodes.(AgentNodeName).fluxes./sum(Nodes.ListNodes.(AgentNodeName).fluxes);
    
    %Randomly select one of the children based on their probability
    Agents.(char(ListOfAgents(i))).currentNode = char(Nodes.ListNodes.(AgentNodeName).children(find(rand<cumsum(Nodes.ListNodes.(AgentNodeName).probabilities),1,'first')));
end

