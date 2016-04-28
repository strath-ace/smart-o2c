function [newnode_ID,childID] = ChooseNode(currentNode,posschildren)
% This function is used to choose a new node for ramification
%
% Inputs:
% * posschildren   : Vector containing the node_IDs of all the possible nodes
% 
% Outputs:
% * newnode_ID        : The node_ID chosen
%
% Author: Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk

temp = strsplit(currentNode,'____');
if length(temp)==1
    exchild = temp{1};
else
    exchild = temp{2};
end

%Choose a random decision (node_ID)
randdecision = randi([1 length(posschildren)]);
childID = posschildren{randdecision};
newnode_ID = strcat(exchild,'____',childID);


end

