function [newNode] = CreateNode(Inputs,ListNodes,node_ID,parent)
%This function creates the structure for a new node.
%
% Inputs:
% * Inputs         : Structure containing the PhysarumSolver inputs
% * ListNodes      : Structure containing the graph
% * decisionname   : string with the chosen target
% * characteristic : a number that describes this new node (eg ToF)
% * parent         : a string containing the node_ID of the parent node
%
% Outputs: 
% * newNode : a structure describing the new node and its attributes
%
% Author: Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk

        
%Split the newnode_ID into the chosen target & characteristic       
temp = strsplit(node_ID, '___');
decisionname = char(temp(1)); 

charstring = strsplit(char(temp(2)),'__');
charwithdot = strrep(charstring,'_','.');

characteristics = str2double(charwithdot);

%Find the decision that was made by the parent
parentdecision = strsplit(parent,'_');
parentdecision = parentdecision(1);

%Add the parent's decision to the list of previous decisions
previousdecisions = [ListNodes.(parent).previousdecisions, {parentdecision}]; 
previousdecisions = [previousdecisions{:}];

%Determine the possible decisions & the number of times each target can
%still be visisted
[possibledecisions, visitsleft] = DeterminePossDecisions(Inputs, ListNodes, parent, previousdecisions, decisionname);

%Create structure of the new node
newNode = struct('node_ID',           node_ID,... % The ID of the node
                 'parent',            parent, ... % The parent of the node
                 'children',          [],... % Matrix that holds the nodes' connections to each other
                 'radius',            [Inputs.StartingRadius],... % The radius of each connection
                 'length',            [],... % The length of each connection
                 'flux',              [],... % Matrix containing each connection's flux
                 'characteristics',   [SetNodeAttributes(Inputs,characteristics)], ... % Characteristics that describe this node (such as orbital elements & ToF .)
                 'previousdecisions', {previousdecisions},... %List of previous decisions made                    
                 'possibledecisions', {possibledecisions}, ... %List containing the decisions that can still be made
                 'VisitsLeft',        {visitsleft} ... % Vector containing the number of times each target cna still be visisted
                );
            
 
%Add the length of the structure. This can only be done after the creation of the structure, as the CostFunction itself needs it             
newNode.length = Inputs.CostFunction(ListNodes.(parent), newNode);

%prevent the length from being 0 (and the flux from becoming inf)
if (newNode.length == 0)
    newNode.length = Inputs.IfZeroLength;
end 

%Calculate the flux and add it to the structure
newNode.flux = [CalculateFlux(Inputs,newNode)];
                 
end
