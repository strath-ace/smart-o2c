function [Agents] = CreateAgents(ListNodes, NumberOfAgents)
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%-----------Copyright (C) 2018 University of Strathclyde and Authors-----------
%
%
%
%% CreateAgents: This function creates the structure containing the Agent decisions and other characteristics.
% 
%% Inputs:
% * ListNodes      : Structure containing the graph
% * NumberOfAgents : The number of agents used in the solver [integer]
% 
%% Outputs: 
% * Agents         : The structure containing the set of agents and their
%                    characteristics
% 
%% Author(s): Aram Vroom (2016)
% Email:  aram.vroom@strath.ac.uk

% Initialise the Agents structure
Agents = struct();          

% Retreive the currently existing nodes and the name of the root
existingnodes = fieldnames(ListNodes);
rootname = existingnodes(1);

% Loop over the number of agents to be created
for i = 1:NumberOfAgents

    % Find the field name of the specific agent in the Agents structure
    agentid = ['agent', num2str(i)];

    % Create the structure for each agent
    AgentStruct = struct('previousListNodes',       [], ...
                         'previouscosts',           [], ...
                         'currentNode',       [rootname] ...
                        );

    % Add agent's structure to the Agents main structure                              
    Agents.(agentid) = AgentStruct;                            

end

