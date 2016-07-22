function [ListNodes] = Dilation(Inputs, ListNodes, Agents, agent)
%% Dilation: This function handles the dilation and evaporation of the paths taken by an agent
%
%% Inputs:
% * Inputs      : Structure containing the PhysarumSolver inputs
% * ListNodes   : Structure containing the graph
% * Agents      : Structure containing the Agents
% * agent       : Cell with the current agents' name
% 
%% Outputs:
% * ListNodes   : Structure containing the updated graph
%
%% Author: Aram Vroom (2016)
% Email:  aram.vroom@strath.ac.uk

%Convert the agent input into a character array
agent = char(agent);

%List all the nodes the agent has visisted
visistednodes = [Agents.(agent).previousListNodes {Agents.(agent).currentNode}];

%Calculate the total cost of the path taken by the agent
totalcost = sum(Agents.(agent).previouscosts);

%Loop over each of the nodes visisted
for i = 1:length(visistednodes)
    
    %For ease of reading, define the node currently being evaluated & its parent
    %as a separate variable
    evaluatednode = char(visistednodes(i));
    parent = ListNodes.(evaluatednode).parent;
    
    %Check if the node has a parent - ignores root
    if ~isempty(parent)
        
        %Dilate the respective link in the parent's node structure
        ListNodes.(evaluatednode).radius = ListNodes.(evaluatednode).radius + Inputs.LinearDilationCoefficient*ListNodes.(evaluatednode).radius/totalcost;
        
         %Check if the link's radius is not too large or small. Correct if
        %so, set to max/min radius
        ListNodes.(evaluatednode).radius(ListNodes.(evaluatednode).radius./Inputs.StartingRadius > Inputs.MaximumRadiusRatio) = Inputs.MaximumRadiusRatio*Inputs.StartingRadius;
        ListNodes.(evaluatednode).radius(ListNodes.(evaluatednode).radius./Inputs.StartingRadius < Inputs.MinimumRadiusRatio) = Inputs.MinimumRadiusRatio*Inputs.StartingRadius;
        
        %Update flux
        ListNodes.(evaluatednode).flux = CalculateFlux(Inputs,ListNodes.(evaluatednode));
    end

end

end

       
