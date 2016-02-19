function [Nodes] = AddNode(Nodes)
% This function adds a node to the existing graph
%   
% Inputs:
% * Nodes: Structure that contains the currently existing nodes
%
% OUTPUTS: 
% * Nodes: Previous Nodes structure but with additional node
%
% Author: Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk

%Find the current size of the matrix
currentsize = size(Nodes.links);

%Define the the size it becomes when adding a row and column
newsize = currentsize+1;

%Obtain the fieldnames of the Nodes structure
fields = fieldnames(Nodes);

%Loop over all the fields of structure
for i = 1:length(fields);
    
    %Define an empty matrix of the new size
    newmatrix = zeros(newsize);

    %Paste the current matrix of the respective field in the Node structure
    %in the top left corner of this empty matrix
    newmatrix(1:currentsize,1:currentsize)=Nodes.(fields{i});
    
    %Set the matrix in the Node structure equal to this matrix with the
    %additional row and column
    Nodes.(fields{i}) = newmatrix;
end
 