function [Solutions, BestSolution, InitializedInputs, ListNodes, Agents] = PhysarumSolver(InitializedInputs, ListNodes)
% This script contains the main logic of AIDMAP solver. 
%
% Inputs:
% * InitializedInputs  : The structure containing the options set by the
%                        user
% * ListNodes          : Structure containing the initial list of nodes
%
% Outputs: 
% * Solutions          : The structure containing the solutions found
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

BestSolution.BestChain = [];
BestSolution.BestCost = [];

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
        
        %Update veins with the dilation and evaporation mechanics
        [ListNodes] = DilationEvaporation(InitializedInputs, ListNodes, Agents, agentnames(i));
        
    %End agent loop
    end
    
    %Update the veins with the growth factor mechanic
    [ListNodes, BestSolution] = GrowthFactor(InitializedInputs, ListNodes, Solutions, BestSolution);
    
    %Check whether the algorithm should be restarted
    restartflag = RestartCheck(InitializedInputs, Agents);
    
    %If so, reset the veins 
    if (restartflag && j ~= InitializedInputs.Generations) %~= as check whether it works (doesn't reset when last generation is completed)
        ListNodes = RadiusFluxReset(InitializedInputs, ListNodes);
    end

%End generation loop
end


end

