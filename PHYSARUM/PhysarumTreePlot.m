function [] = PhysarumTreePlot(ListNodes)
%This function plots the tree corresponding to the physarum graph
%
% Inputs:
% * ListNodes       : Structure containing the graph
%
% Outputs: 
% * 
%
% Author: Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk

%Retireve the node names
nodenames = fieldnames(ListNodes);

%Loop over all the nodes
for i = 2:length(nodenames)
    %Find each node's parent
    treeplotnodes(i) = strmatch(char(ListNodes.(char(nodenames(i))).parent), nodenames, 'exact');
end

%Plot the tree
treeplot(treeplotnodes);

%Find the x & y locations of each node
[x,y] = treelayout(treeplotnodes);
x = x';
y = y';

%Set the node_ID as label for each node
text(x(:,1), y(:,1), nodenames, 'FontSize', 16, 'VerticalAlignment', 'bottom',...
    'HorizontalAlignment', 'right', 'Interpreter', 'none');

%Remove labels & axis ticks
set(gca, 'XTick', [])
set(gca, 'YTick', [])
xlabel('')
ylabel('')


end

