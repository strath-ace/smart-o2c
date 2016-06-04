function [newnode_ID,nodeindex] = ChooseNode(currentNode,posschildren)
% This function is used to choose a new node for ramification
%
% Inputs:
% * currentNode    : The unique ID of the current node
% * posschildren   : Vector containing the node_IDs of all the possible nodes
% 
% Outputs:
% * newnode_ID     : The node_ID chosen
% * nodeindex      : The index of the node in the posschildren array
%
% Author: Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk

%Split the ID of the curren node
temp = strsplit(currentNode,'____');

%Find the current target + ID attributes. If the current node is the root, 
%the length of the temp vector is 1, and the current node is therefore the root. 
%If this is not the root, the current target is second part of the vector.
if length(temp)==1
    exchild = temp{1};
else
    exchild = temp{2};
end

%Choose a random decision (node_ID)
nodeindex = randi([1 length(posschildren)]);
childID = posschildren{nodeindex};
newnode_ID = strcat(exchild,'____',childID);


end

