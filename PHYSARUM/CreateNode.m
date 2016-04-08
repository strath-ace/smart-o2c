function [newNode] = CreateNode(Inputs,Nodes,decisionname,characteristic,parent)
%This function creates the structure for a new node.
%
% Inputs:
% * Inputs         : Structure containing the PhysarumSolver inputs
% * Nodes          : Structure containing the graph
% * decisionname   : string with the chosen target
% * characteristic : a number that describes this new node (eg ToF)
% * parent         : a string containing the node_ID of the parent node
%
% Outputs: 
% * newNode : a structure describing the new node and its attributes
%
% Author: Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk

%Generate the node ID using the decision name (eg target) & the node's
%characteristic
node_ID = strcat(decisionname,{'_'},num2str(characteristic));

%Find the decision that was made by the parent
parentdecision = strsplit(parent,'_');
parentdecision = parentdecision(1);

%Add the parent's decision to the list of previous decisions
previousdecisions = [Nodes.(parent).previousdecisions, {parentdecision}]; 
previousdecisions = [previousdecisions{:}];

%Determine the possible decisions & the number of times each target can
%still be visisted
[possibledecisions,visitsleft] = DeterminePossDecisions(Inputs,Nodes,parent,previousdecisions,decisionname);

%Create structure of the new node
newNode = struct('node_ID',           node_ID,... % The ID of the node
                 'parent',            parent, ... % The parent of the node
                 'children',          [],... % Matrix that holds the nodes' connections to each other
                 'radius',            [],... % The radius of each connection
                 'pressure_gradient', [],... % The pressure gradient over each connection
                 'lengths',           [],... % The length of each connection
                 'fluxes',            [],... % Matrix containing each connection's flux
                 'probabilities',     [],... % Matrix containing the probability for each connection
                 'characteristics',   [characteristic], ... % Characteristics that describe this node (such as orbital elements & ToF .)
                 'previousdecisions', {previousdecisions},... %List of previous decisions made                    
                 'possibledecisions', {possibledecisions}, ... %List containing the decisions that can still be made
                 'VisitsLeft',        {visitsleft} ... % Vector containing the number of times each target cna still be visisted
                );
             

                 
end

