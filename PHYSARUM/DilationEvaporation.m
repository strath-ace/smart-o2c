function [Nodes] = DilationEvaporation(Inputs, Nodes, Agents, agent)
%This function handles the dilation and evaporation of the paths taken by
%an agent
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

%Convert the agent input into a character array
agent = char(agent);

%List all the nodes the agent has visisted
visistednodes = [Agents.(agent).previousNodes {Agents.(agent).currentNode}];

%Calculate the total cost of the path taken by the agent
totalcost = sum(Agents.(agent).previouscosts);

%Loop over each of the nodes visisted
for i = 1:length(visistednodes)
    
    %For ease of reading, define the node currently being evaluated & its parent
    %as a separate variable
    evaluatednode = char(visistednodes(i));
    parent = Nodes.(evaluatednode).parent;
    
    %Check if the node has a parent
    if ~isempty(parent)
        
        %Find the node's index in the children list of the parent
        indexofnode = strcmp(Nodes.(parent).children,evaluatednode);
        
        %Dilate the respective link in the parent's node structure
        Nodes.(parent).radius = Nodes.(parent).radius + Inputs.LinearDilationCoefficient*indexofnode.*Nodes.(parent).radius/totalcost;
        
        %Simulate evaporation in this link
        Nodes.(parent).radius = Nodes.(parent).radius - Inputs.EvaporationCoefficient*indexofnode.*Nodes.(parent).radius;

        %Check if the link's radius is not too large or small. Correct if
        %so
        Nodes.(parent).radius(Nodes.(parent).radius > Inputs.MaximumRadius) = Inputs.MaximumRadius;
        Nodes.(parent).radius(Nodes.(parent).radius < Inputs.MinimumRadius) = Inputs.MinimumRadius;
    end

end

