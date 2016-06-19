function [Solutions, BestSolution, InitializedInputs, ListNodes, Agents, History] = PhysarumSolver(InitializedInputs, ListNodes)
% This script contains the main logic of AIDMAP solver. 
%
% Inputs:
% * InitializedInputs  : The structure containing the options set by the
%                        user
% * ListNodes          : Structure containing the initial list of nodes
%
% Outputs: 
% * Solutions          : The structure containing the solutions found
% * BestSolution       : The best solution found
% * InitializedInputs  : The structure containing the options set by the
%                        user
% * ListNodes          : Structure containing the final structure with the
%                        nodes
% * Agents             : the structure containing the set of agents and their
%                        characteristics
%
% Author: Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk

%Initialize the Solutions structure
Solutions.Nodes = [];
Solutions.Costs = [];

%Initialize the structure that will contain the best solution
BestSolution.BestChain = [];
BestSolution.BestCost = [];
History.radius = {};
History.BestSolution = {};
History.BestCost = {};
History.AgentMovement = {};

%Loop over the generations
for j = 1:InitializedInputs.Generations
    
    %Show the current generation
    disp(strcat(['Starting Generation',' ',num2str(j)]));
    
    %Create a new structure for the agents
    Agents = CreateAgents(ListNodes,InitializedInputs.NumberOfAgents);
    
    %Retrieve the agent names
    agentnames = fieldnames(Agents);
    
    %Loop over the agents
    for i = 1:length(agentnames);
        
        %Show the current agent
        disp(strcat(['Moving Agent',' ',num2str(i)]));
        
        %Reset the agent death flag
        agentdeathflag = 0;
        
        %Continue moving the agent until the death flag becomes 1
        while ~agentdeathflag
            [Solutions, ListNodes, Agents, agentdeathflag] = AgentMovement2(InitializedInputs, Solutions, ListNodes, Agents, agentnames(i));
                  
            
        end
        if ((InitializedInputs.SaveHistory ~= 0) ||(InitializedInputs.GenerateGraphPlot ~= 0))
        nodenames = fieldnames(ListNodes);     
        for p = 2:length(nodenames)
                radii(p) = ListNodes.(char(nodenames(p))).radius;
        end
        History.radius(end+1) = {radii};
        History.AgentMovement{j,i} = [Agents.(char(agentnames(i))).previousListNodes Agents.(char(agentnames(i))).currentNode];
        end
        
        %Update veins with the dilation and evaporation mechanics
        [ListNodes] = Dilation(InitializedInputs, ListNodes, Agents, agentnames(i));
    
                       
    %End agent loop
    end
    
   
    %Update the veins with the growth factor mechanic
    [ListNodes, BestSolution] = GrowthEvaporation(InitializedInputs, ListNodes, Solutions, BestSolution);
    
    if ((InitializedInputs.SaveHistory ~= 0) ||(InitializedInputs.GenerateGraphPlot ~= 0))
        History.BestSolution(end+1) = BestSolution.BestChain(1);
        History.BestCost(end+1) = BestSolution.BestCost(1);
    end
    
    %Check whether the algorithm should be restarted
    restartflag = RestartCheck(InitializedInputs, Agents);
    
    %If so, reset the veins 
    if (restartflag)
        ListNodes = RadiusFluxReset(InitializedInputs, ListNodes);
    end

%End generation loop
end


end

