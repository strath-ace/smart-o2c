function [Agents] = CreateAgents(ListNodes, NumberOfAgents)
%% CreateAgents: This function creates the structure containing the Agent decisions and other characteristics.
%
%% Inputs:
% * ListNodes      : Structure containing the graph
% * NumberOfAgents : The number of agents used in the solver
%
%% Outputs: 
% * Agents         : The structure containing the set of agents and their
%                    characteristics
%
%% Author: Aram Vroom (2016)
% Email:  aram.vroom@strath.ac.uk

%Initialise the Agents structure
Agents = struct();          

%Retreive the currently existing nodes and the name of the root
existingnodes = fieldnames(ListNodes);
rootname = existingnodes(1);

%Loop over the number of agents to be created
for i = 1:NumberOfAgents

    %Find the field name of the specific agent in the Agents structure
    agentid = ['agent', num2str(i)];

    %Create the structure for each agent
    AgentStruct = struct('previousListNodes',       [], ...
                         'previouscosts',           [], ...
                         'currentNode',       [rootname] ...
                        );

    %Add agent's structure to the Agents main structure                              
    Agents.(agentid) = AgentStruct;                            

end

