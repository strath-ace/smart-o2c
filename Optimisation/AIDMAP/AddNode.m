function [ListNodes] = AddNode(ListNodes, newNode)
% This function adds a node to the list of nodes
%   
% Inputs:
% * ListNodes    : Structure that contains the currently existing nodes
% * newNode  : Structure that cotnains the new node and its attributes
%
% OUTPUTS: 
% * ListNodes: Previous ListNodes structure but with additional node
%
% Author: Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk

%Retrieve the node_ID of the node that is added
newnode_ID = newNode.node_ID;

%Add the new node underneath the currently existing stucture
ListNodes.(newnode_ID) = newNode;

%Retrieve the node_ID of the node's parent
parent = ListNodes.(newnode_ID).parent;

%Check if the noe has a parent
if ~isempty(parent)
    
    %Add the node to the children list of the parent node in the the ListNodes 
    %structure
    ListNodes.(parent).children = [ListNodes.(parent).children {newnode_ID}];  

end
    

end