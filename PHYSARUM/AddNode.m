function [Nodes] = AddNode(Inputs,Nodes,newNode,linkcost)
% This function adds a node to the existing graph
%   
% Inputs:
% * Inputs   : Structure containing the PhysarumSolver inputs
% * Nodes    : Structure that contains the currently existing nodes
% * newNode  : Structure that cotnains the new node and its attributes
% * linkcost : The cost of the connection between the parent and the new
%              node
%
% OUTPUTS: 
% * Nodes: Previous Nodes structure but with additional node
%
% Author: Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk

%Retrieve the node_ID of the node that is added
newnode_ID = newNode.node_ID;

%Add the new node underneath the currently existing stucture
Nodes.(newnode_ID) = newNode;

%Retrieve the node_ID of the node's parent
parent = Nodes.(newnode_ID).parent;

%Check if the noe has a parent
if ~isempty(parent)
    
    %Add the node to the children list of the parent node in the the Nodes 
    %structure
    Nodes.(parent).children = [Nodes.(parent).children {newnode_ID}];  
    
    %Add the link's cost to the parent node
    Nodes.(parent).lengths = [Nodes.(parent).lengths linkcost];
    Nodes.(parent).radius = [Nodes.(parent).radius Inputs.StartingRadius];
end
    

end