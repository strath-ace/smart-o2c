function [Nodes,Agents] = Translation(Nodes,Agents,agent)
% This function is used to move agents through the graph/tree
%
% Inputs:
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

currentagent = char(agent);
     
%For easy of reading, define the name of the agent's current node
%as a separate variable
agentnodename =  char(Agents.(currentagent).currentNode);

%Add agent's current node to its history of previous nodes
Agents.(currentagent).previousNodes = [Agents.(currentagent).previousNodes {Agents.(currentagent).currentNode}];

%Calculate the probabilities for the children using their fluxes
Nodes.ListNodes.(agentnodename).probabilities = Nodes.ListNodes.(agentnodename).fluxes./sum(Nodes.ListNodes.(agentnodename).fluxes);

%Randomly select one of the children based on their probability
Agents.(currentagent).currentNode = char(Nodes.ListNodes.(agentnodename).children(find(rand<cumsum(Nodes.ListNodes.(agentnodename).probabilities),1,'first')));


