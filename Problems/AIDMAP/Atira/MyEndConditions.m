function [continueflag] = MyEndConditions(Inputs, Agents, agent)
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
%% MyEndConditions: This function confirms whether the final condtions have been reached
% 
%% Inputs:
% * Inputs        : Structure containing the PhysarumSolver inputs
% * Agents        : Structure containing the agents
% * agent         : The name of the agent currently being evaluated [string]
% 
%% Outputs: 
% * continueflag  : The flag that confirms whether the algorithm has
%                  reached the final condition. This should be 1 if 
%                  it has,  and 0 otherwise
% 
%% Author(s): Aram Vroom (2016)
% Email:  aram.vroom@strath.ac.uk

% Retrieve the current node
currentnode = Agents.(char(agent)).currentNode;

% Obtain the final conditions
endconditions = Inputs.EndConditions;

% Check if final city reached
temp = strsplit(currentnode, '_');
currentcity = temp{end-3};
check1 = ~(sum(ismember(currentcity,  endconditions{1}))==length(currentcity));

% Determine the continueflag
continueflag = check1;

end

