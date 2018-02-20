function [Solutions, BestSolution, ListNodes, Agents, History, funccalls] = PhysarumSolver(InitialisedInputs, ListNodes)
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
<<<<<<< HEAD
%-----------Copyright (C) 2018 University of Strathclyde and Authors-----------
=======
%-----------Copyright (C) 2016 University of Strathclyde-------------
>>>>>>> 5b7361d93c9119cf1d2e9e6c885bed93f924d71b
%
%
%
%% PhysarumSolver: This script contains the main logic of the AIDMAP solver. 
% 
%% Inputs:
% * InitialisedInputs  : The structure containing the options set by the user
% * ListNodes          : Structure containing the initial list of nodes
% 
%% Outputs: 
% * Solutions          : The structure containing the solutions found
%                      * Solutions.Nodes: cell array containing all the
%                        solutions (paths) found
%                      * Solutions.Costs: cell array containing the costs
%                        corresponding to each link of each solution found
% * BestSolution       : Structure containing the best solution found
%                      * BestSolution.BestChain = cell array containing
%                        the node IDs of the best chain
%                      * BestSolution.BestCost = the total cost
%                        corresponding to the best chain               
% * ListNodes          : Structure containing the final structure with the
%                        nodes
% * Agents             : The structure containing the set of agents and their
%                        characteristics
% * History            : The vein radii and agent movement throughout the
%                        generations [structure]
% * funccalls          : The number of cost function calls [integer]
% 
%% Author(s): Aram Vroom (2016)
% Email:  aram.vroom@strath.ac.uk

% Initialise the Solutions structure
Solutions.Nodes = [];
Solutions.Costs = [];

% Set the initial number of function calls
funccalls = 0;

% Initialise the structure that will contain the best solution
BestSolution.BestChain = [];
BestSolution.BestCost = [];

% Initialzie the structure that will contain the history
History.radius = {};
History.BestSolution = {};
History.BestCost = {};
History.AgentMovement = {};

% Loop over the generations
for j = 1:InitialisedInputs.Generations
    
    % Show the current generation
    disp(strcat([datestr(now), ' === Starting Generation', ' ', num2str(j)]));
    
    % Create a new structure for the agents
    Agents = CreateAgents(ListNodes, InitialisedInputs.NumberOfAgents);
    
    % Retrieve the agent names
    agentnames = fieldnames(Agents);
    
    % Loop over the agents
    for i = 1:length(agentnames);
        
        % Show the current agent
        disp(strcat([datestr(now), ' === Moving Agent', ' ', num2str(i)]));
        
        % Reset the agent death flag
        agentdeathflag = 0;
        
        % Continue moving the agent until the death flag becomes 1
        while ~agentdeathflag
            [Solutions, ListNodes, Agents, agentdeathflag, funccalls] = AgentMovement(InitialisedInputs, Solutions, ListNodes, Agents, agentnames(i), funccalls);
                  
            
        end        
        
        % Update veins with the dilation and evaporation mechanics
        [ListNodes] = Dilation(InitialisedInputs, ListNodes, Agents, agentnames(i));
        
        % If the user has specified the request to save the history or to
        % generate the graph plot, save the current radii of the veins and
        % the path that each agent has moved
        if ((InitialisedInputs.SaveHistory ~= 0) ||(InitialisedInputs.GenerateGraphPlot ~= 0))
        
            % Obtain the names of the currently existing nodes
            nodenames = fieldnames(ListNodes); 
            
            % Check if additional nodes apart from the root have been found
            if length(nodenames)> 1   
                
                % If so, loop over all the nodes and save their radius
                for p = 2:length(nodenames)
                        radii(p) = ListNodes.(char(nodenames(p))).radius;
                end
                
                % Save the radius and the agent movement
                History.radius(end+1) = {radii};
                History.AgentMovement{j, i} = [Agents.(char(agentnames(i))).previousListNodes Agents.(char(agentnames(i))).currentNode];
            end
        end
                          
    end
    
   
    % Update the veins with the growth factor mechanic
    [ListNodes, BestSolution] = GrowthEvaporation(InitialisedInputs, ListNodes, Solutions, BestSolution);
    
    % If the user has specified to save the history or the generate the
    % graph lot, save the best solution and best cost
    if ((InitialisedInputs.SaveHistory ~= 0) ||(InitialisedInputs.GenerateGraphPlot ~= 0))
        History.BestSolution(end+1) = BestSolution.BestChain;
        History.BestCost(end+1) = BestSolution.BestCost;
    end
    
    % Check whether the algorithm should be restarted
    restartflag = RestartCheck(InitialisedInputs, Agents);
    
    % If so, reset the veins 
    if (restartflag)
        ListNodes = RadiusFluxReset(InitialisedInputs, ListNodes);
        
        % If the user has specified the request to save the history or to
        % generate the graph plot, save the current radii of the veins and
        % the path that each agent has moved
        if ((InitialisedInputs.SaveHistory ~= 0) ||(InitialisedInputs.GenerateGraphPlot ~= 0))
            
            % Obtain the names of the currently existing nodes
            nodenames = fieldnames(ListNodes); 
            
            % Check if additional nodes apart from the root have been found
            if length(nodenames)> 1
                
                % Loop over all the nodes
                for p = 2:length(nodenames)
                    
                    % Save the radius
                    radii(p) = ListNodes.(char(nodenames(p))).radius;
                end
                
                % Add the radii and the agent movement to their respective cells
                History.radius(end) = {radii};
                History.AgentMovement{j, i} = [Agents.(char(agentnames(i))).previousListNodes Agents.(char(agentnames(i))).currentNode];
            end
        end
    end

end


end

