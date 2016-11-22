function [ListNodes] = CreateListNodes(Inputs)
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%-----------Copyright (C) 2016 University of Strathclyde-------------
%
%
%
%% CreateListNodes: This function creates the structure that holds all information of every created node
% 
%% Inputs:
% * Inputs : Structure containing the PhysarumSolver inputs
% 
%% Outputs: 
% * ListNodes  : Empty structure that will contain the Graph's information
% 
%% Author(s): Aram Vroom (2016)
% Email:  aram.vroom@strath.ac.uk


%% Create the Unique Identifier (UID)
% As only underscores can be used in field names, the UID is made up as follows:
% 2 underscores define the difference between the city and the attributes
% 1 underscore defines the difference between two attributes.

rootattribstr = [];
for i = 1:(length(Inputs.RootAttrib)*2+1)
    if rem(i, 2)
        rootattribstr = strcat(rootattribstr, '0');
    else
        rootattribstr = strcat(rootattribstr, '_');
    end
end

rootID = strcat(Inputs.RootName, '__', rootattribstr);

% Create the ListNodes structure
ListNodes = struct(struct(rootID,   ...
                         struct(...
                         'node_ID',             rootID, ...      % The node_ID of the root
                         'parent',              [], ...          % The parent of the node
                         'children',            [], ...          % Matrix that holds the nodes' connections to each other
                         'radius',              [], ...          % The radius of the connection to the parent
                         'length',              [], ...          % The length of the connection to the parent
                         'flux',                [], ...          % The flux of the connection to the parent
                         'attributes',          [SetNodeAttributes(Inputs, [], ...          % Problem-specific Attributes that 
                                                 Inputs.RootName, ...                       % describe this node 
                                                 zeros(1, length(Inputs.RootAttrib)))], ... % (such as orbital elements & Time of Flight. etc)
                         'previousdecisions',    [], ...         % List of the previous decisions made
                         'possibledecisions',    {Inputs.PossibleDecisions}, ...            % Cities that can still be visisted by the node
                         'VisitsLeft',           {Inputs.MaxVisits} ...                     % Vector containing the number of times each city can still be visisted                         
                     )));

% This will track the validity of attempted children                       
if Inputs.LowMem == 0
    ListNodes.(rootID).ChildValidityTracker = 1:length(Inputs.PossibleListNodes);
end

end 
