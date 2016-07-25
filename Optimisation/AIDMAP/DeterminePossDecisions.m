% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%-----------Copyright (C) 2016 University of Strathclyde-------------
%


function [possibledecisions, visitsleft] = DeterminePossDecisions(Inputs, ListNodes, parent, previousdecisions, decisionname, decisionindex)
%% DeterminePossDecisions: This function determines the possible cities the agent can move to from the current city, using the maximum number of (consecutive) visits
% 
%% Inputs:
% * Inputs            : Structure containing the PhysarumSolver inputs
% * ListNodes         : Structure that contains the Graph's information
% * parent            : A string containing the node_ID of the parent node
% * previousdecisions : A list of the previous decisions made
% * decisionname      : String with the name of the chosen city
% * decisionindex     : The index of the city's name in the list the Inputs.Cities array
% 
%% Outputs: 
% * possibledecisions : A list of the decisions that can still be made by
%                       the node
% * visitsleft        : Updated vector with the number of times each city
%                       can still be visisted
% 
%% Author(s): Aram Vroom (2016)
% Email:  aram.vroom@strath.ac.uk


% Create a concatenated list of the previous decisions, including the newly generated node 
previousdecisions = strcat(previousdecisions{:}, decisionname);

% Retrieve the maximum number of consecutive resonance orbits
maxconsecutive = Inputs.MaxConsecutiveVis(decisionindex);

% Copy the number of visits left of the parent
visitsleft = ListNodes.(parent).VisitsLeft;

% Find the number of times a city can still be visited
visitsleft(decisionindex) = ListNodes.(parent).VisitsLeft(decisionindex)-1;

% Retrieve all the possible decisions
alldecisions = Inputs.PossibleDecisions;

% Initially, set the decisions the node can make equal to be equal to all the decisions
possibledecisions = alldecisions;

% Remove decisions with visitsleft = 0
possibledecisions(visitsleft == 0) = [];

% Check whether resonance orbit is allowed by first generating the string that
% should be checked for in the previous decisions
checkforconcities = repmat(decisionname, 1, maxconsecutive+1);

% See if max. resonant orbit has been performed as last trajectory so far 
resperformed = max(strfind(previousdecisions, checkforconcities))==(length(previousdecisions)-length(checkforconcities)+1);

% If so, exclude the respective city from the list of possible decisions
if resperformed 
    decisionindex = strcmp(decisionname, possibledecisions);
    possibledecisions(decisionindex) = [];
end



end

