function [ListNodes] = AddNode(Inputs, ListNodes, newNode, linkcost)
% This function adds a node to the list of nodes
%   
% Inputs:
% * Inputs   : Structure containing the PhysarumSolver inputs
% * ListNodes    : Structure that contains the currently existing nodes
% * newNode  : Structure that cotnains the new node and its attributes
% * linkcost : The cost of the connection between the parent and the new
%              node
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
    
    %Add the link's cost to the parent node
    ListNodes.(parent).lengths = [ListNodes.(parent).lengths linkcost];
    ListNodes.(parent).radius = [ListNodes.(parent).radius Inputs.StartingRadius];
    
    %Add the link's flux to the parent node
    ListNodes.(parent).fluxes = [ListNodes.(parent).fluxes CalculateFlux(Inputs,ListNodes,newnode_ID)];
end
    

end