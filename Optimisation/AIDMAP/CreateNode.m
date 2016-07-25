function [newNode] = CreateNode(Inputs, ListNodes, node_ID, parent)
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%-----------Copyright (C) 2016 University of Strathclyde-------------
%
%
%
%% CreateNode: This function creates the structure for a new node.
% 
%% Inputs:
% * Inputs         : Structure containing the PhysarumSolver inputs
% * ListNodes      : Structure containing the graph
% * node_ID        : The unique identifier of the node
% * parent         : A string containing the node_ID of the parent node
% 
%% Outputs: 
% * newNode : A structure describing the new node and its attributes
% 
%% Author(s): Aram Vroom (2016)
% Email:  aram.vroom@strath.ac.uk
        
% Retrieve the city name and attributes from the unique ID
temp = strsplit(char(node_ID), '_');
cityname = temp{end-length(Inputs.AttributeIDIndex)-1};
attributes = str2double(temp(end-length(Inputs.AttributeIDIndex):end));

% Calculate the attributes and use the cost function to obtain the vein's length
Attributes = SetNodeAttributes(Inputs, ListNodes.(parent), cityname, attributes);
[Attributes, veinlength] = Inputs.CostFunction(Inputs, ListNodes.(parent), Attributes);

% Sanity check
if (veinlength == Inf)
    newNode = [];
    return
end

% If the node is valid according to the first check, use the created
% node structure to further determine the validity
[checktot] = Inputs.CreatedNodeCheckFile(Inputs, Attributes, ListNodes, parent);

% Sanity check #2
if (checktot == 0)
    newNode = [];
    return
end

% Prevent the length from being 0 (and the flux from becoming inf)
veinlength(veinlength == 0) = Inputs.IfZeroLength;

% Find the decision that was made by the parent
parentdecision = strsplit(parent, '_');
parentdecision = parentdecision(1);

% Add the parent's decision to the list of previous decisions
previousdecisions = [ListNodes.(parent).previousdecisions, {parentdecision}]; 
previousdecisions = [previousdecisions{:}];

% Determine the possible decisions & the number of times each city can
% still be visisted
[possibledecisions, visitsleft] = DeterminePossDecisions(Inputs, ListNodes, parent, previousdecisions, cityname, attributes(1));

% Create structure of the new node
newNode = struct('node_ID',           node_ID, ...                   % The ID of the node
                 'parent',            parent, ...                    % The parent of the node
                 'children',          [], ...                        % Matrix that holds the nodes' connections to each other
                 'radius',            [Inputs.StartingRadius], ...   % The radius of each connection
                 'length',            [veinlength], ...              % The length of each connection
                 'flux',              [], ...                        % Matrix containing each connection's flux
                 'attributes',        [Attributes], ...              % Attributes that describe this node (such as orbital elements & Time of Flight etc.)
                 'previousdecisions', {previousdecisions}, ...       % List of previous decisions made                    
                 'possibledecisions', {possibledecisions}, ...       % List containing the decisions that can still be made
                 'VisitsLeft',        {visitsleft} ...               % Vector containing the number of times each city cna still be visisted                 
                ); 

% This will track the validity of attempted children
if Inputs.LowMem == 0
    newNode.ChildValidityTracker = 1:length(Inputs.PossibleListNodes);
end
 
% Calculate the flux and add it to the structure
newNode.flux = [CalculateFlux(Inputs, newNode)];
                 
end


