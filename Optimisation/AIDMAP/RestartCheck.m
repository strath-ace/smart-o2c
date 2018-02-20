function [restartflag] = RestartCheck(Inputs, Agents)
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%-----------Copyright (C) 2018 University of Strathclyde and Authors-----------
%
%
%
%% RestartCheck: This function checks the restarting of the Physarum algorithm
% 
%% Inputs:
% * Inputs   : Structure containing the PhysarumSolver inputs
% * Agents   : Structure containing the Agents & their characteristics
% 
%% Outputs: 
% * restartflag : A flag that is set to 1 if the algorithm is to be
%                 restarted
% 
%% Author(s): Aram Vroom (2016)
% Email:  aram.vroom@strath.ac.uk

% Initialise the restart flag
restartflag = 0;

% Retrieve the agent names
agents = fieldnames(Agents);

% Loop over all the agents
for i = 1:length(agents)
    
    % Change the agent into characters (needed for structures)
    agent = agents(i);
    agent = agent{1};
    
    % Retrieve the nodes visisted by the agent
    visistednodes{i} = [Agents.(agent).previousListNodes {Agents.(agent).currentNode}];
    
end

% Loop over all the visisted nodes
for i = 1:length(visistednodes)
    for j = 1:length(visistednodes)
        
        % Only compare different rows
        if i ~= j
            
            % Find the number of equal nodes in the sequences
            numberequal = sum(ismember(visistednodes{i}(1, :), visistednodes{j}(1, :)));
            
            % If this exceeds the set threshold, set the restartflag to 1
            if numberequal > Inputs.MinCommonNodesThres
                restartflag = 1;
                break
            end
        end
    end
end

% Show if a restart is needed
if restartflag
    disp('Restart Needed')
end
        
end
