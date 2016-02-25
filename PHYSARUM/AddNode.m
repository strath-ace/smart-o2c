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

%Check whether this is the first node to be added
if isempty(Nodes.ListNodes)
    
    %If so, set the Nodes.ListNodes structure equal to the newNode. The
    %first node is separated due to catstruct not being able to handle
    %empty cells (Nodes.ListNodes)
    Nodes.ListNodes = newNode;
else
           
    %If not, add the new node underneath the currently existing stucture.
    tmp = [fieldnames(Nodes.ListNodes), struct2cell(Nodes.ListNodes); fieldnames(newNode), struct2cell(newNode)].';
    Nodes.ListNodes = struct(tmp{:});
end
 

end