function [ListNodes] = AddNewNode(ListNodes, newNode)
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
%% AddNewNode: This function adds a node to the list of nodes
%  
%% Inputs:
% * ListNodes   : Structure that contains the currently existing nodes
% * newNode     : Structure that cotnains the new node and its attributes
% 
%% Outputs: 
% * ListNodes: Previous ListNodes structure but with additional node
% 
%% Author(s): Aram Vroom (2016)
% Email:  aram.vroom@strath.ac.uk

% Retrieve the node_ID of the node that is added
newnode_ID = newNode.node_ID;

% Add the new node underneath the currently existing stucture
ListNodes.(newnode_ID) = newNode;

% Retrieve the node_ID of the node's parent
parent = ListNodes.(newnode_ID).parent;

% Check if the node has a parent (excludes the root)
if ~isempty(parent)
    
    % If node has a parent, add it to the children list of the parent node 
    % in the the ListNodes structure
    ListNodes.(parent).children = [ListNodes.(parent).children {newnode_ID}];  

end
    

end
