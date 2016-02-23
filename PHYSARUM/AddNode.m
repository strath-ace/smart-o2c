function [Nodes] = AddNode(Inputs,Nodes,parent)
% This function adds a node to the existing graph
%   
% Inputs:
% * Inputs : Structure containing the PhysarumSolver inputs
% * Nodes  : Structure that contains the currently existing nodes
% * parent : The parent of the node to be added
%
% OUTPUTS: 
% * Nodes: Previous Nodes structure but with additional node
%
% Author: Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk

%Find the row to which elements should be added
addrow = length(Nodes.identifier)+1;

%Set the identifier of the added node. Currently equal to the node's row in
%the Nodes.identifier list. Can be changed.
node_ID = addrow;

%Add the (starting) characteristics of each vein to the Nodes structure if
%the Node doesn't have any connections yet
if isempty(Nodes.children{addrow,1})
    Nodes.identifier{addrow,1}        = node_ID;
    
    %Check if this is not the first node.
    if ~isempty(Nodes.parent)
        %Add the parent if this is not the first node
        Nodes.parent{addrow,1}        = parent;

        %Add the generated node to the relevant vectors of the parent
        Nodes.children{parent,1}      = [Nodes.children{parent,1} node_ID];

    else
        %If it is the first node, set the listed parent node_ID to 0
        Nodes.parent{addrow,1}        = 0;         

    end
      
    %Set the initial conditions for the "children" to "probabilities"
    %fields.
    fields = fieldnames(Nodes);
    for i = 3:length(fields)
        Nodes.(fields{i}){addrow,1}       = [];
    end
end

end


