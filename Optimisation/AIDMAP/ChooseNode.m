% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%-----------Copyright (C) 2016 University of Strathclyde-------------
%


function [newnode_ID, nodeindex] = ChooseNode(Inputs, currentNode, posschildren, indextracker, nantracker)
%% ChooseNode: This function is used to choose a new node for ramification
% 
%% Inputs:
% * Inputs                : Structure containing the PhysarumSolver inputs
% * currentNode           : The unique identifier of the current node [string]
% * posschildren          : Vector containing the node_IDs of all the possible nodes
% * indextracker          : The vector containing the data as to whether
%                           each possible child to this node is found to be invalid or not
% * nantracker            : The number of invalid nodes found so far [integer]
% 
%% Outputs:
% * newnode_ID     : The unique identifier of the chosen node [string]
% * nodeindex      : The index of the node in the posschildren array [integer]
% 
%% Author(s): Aram Vroom (2016)
% Email:  aram.vroom@strath.ac.uk

% Split the ID of the curren node
temp = strsplit(currentNode, '___');

% Find the current city + ID attributes. If the current node is the root, 
% the length of the temp vector is 1, and the current node is therefore the root. 
% If this is not the root, the current city is second part of the vector.
if length(temp)==1
    exchild = temp{1};
else
    exchild = temp{2};
end

% Find the total number of possible children
numberofchildren = length(posschildren);

% Calculate the probability of picking a feasible (non-NaN) child
pickprob = 1-(nantracker)/numberofchildren;

% Check whether the probability is larger than the probability specified
if pickprob > Inputs.MinPickProbability
    
    % If the probability is larger, attempt to pick a child by randomly
    % selecting one from the ChildValidityTracker list
    nodeindex = indextracker(randi([1 numberofchildren]));
    
    % If the child node was already found to be infeasible (the
    % ChildValidityTracker is NaN), try again
    while isnan(nodeindex)
        nodeindex = indextracker(randi([1 numberofchildren]));
    end
else
    
    % If the probability is smaller than specified, find index of non-NaN values
    nonanindices = indextracker(~isnan(indextracker));

    % Choose a random decision (node_ID)
    nodeindex = nonanindices(randi([1 length(nonanindices)]));
end

% Retrieve the ID of the child
childID = posschildren{nodeindex};

% Generate the new node ID
newnode_ID = strcat(exchild, '___', childID);
end

