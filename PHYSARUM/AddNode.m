function [Nodes] = AddNode(Nodes,newNode)
% This function adds a node to the existing graph
%   
% Inputs:
% * Nodes   : Structure that contains the currently existing nodes
% * newNode : Structure that cotnains the new node and its attributes
%
% OUTPUTS: 
% * Nodes: Previous Nodes structure but with additional node
%
% Author: Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk

%Add the new node underneath the currently existing stucture
tmp = [fieldnames(Nodes.ListNodes), struct2cell(Nodes.ListNodes); fieldnames(newNode), struct2cell(newNode)].';
Nodes.ListNodes = struct(tmp{:});

%Retrieve the node_ID of the node that is added
node_ID = fieldnames(newNode);

%Retrieve the node_ID of the node's parent
parent = Nodes.ListNodes.(node_ID{1}).parent;
if ~isempty(parent)
        %Add the node to the children list of the parent node in the
        %ListNodes structure within the Nodes structure
        Nodes.ListNodes.(parent).children = [Nodes.ListNodes.(parent).children node_ID];
        
end
    

end