function [newNode] = CreateNode(Inputs, ListNodes, node_ID, parent)
%This function creates the structure for a new node.
%
% Inputs:
% * Inputs         : Structure containing the PhysarumSolver inputs
% * ListNodes      : Structure containing the graph
% * node_ID        : The unique ID of the node
% * parent         : a string containing the node_ID of the parent node
%
% Outputs: 
% * newNode : a structure describing the new node and its attributes
%
% Author: Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk
        
%Retriev ethe target name and attributes from the unique ID
temp = strsplit(char(node_ID),'_');
targetname = temp{end-length(Inputs.AttributeIDIndex)-1};
attributes = str2double(temp(end-length(Inputs.AttributeIDIndex):end));

%Calculate the attributes and use the cost function to obtain the vein's length
Attributes = SetNodeAttributes(Inputs, ListNodes.(parent), targetname, attributes);
[Attributes, veinlength] = Inputs.CostFunction(Inputs, ListNodes.(parent), Attributes);

%Sanity check
if (veinlength == Inf)
    newNode = [];
    return
end

%If the node is valid according to the first check, use the created
%node structure to further determine the validity
[checktot] = Inputs.CreatedNodeCheckFile(Inputs, Attributes, ListNodes, parent);

%Sanity check #2
if (checktot == 0)
    newNode = [];
    return
end

%Prevent the length from being 0 (and the flux from becoming inf)
veinlength(veinlength == 0) = Inputs.IfZeroLength;

%Find the decision that was made by the parent
parentdecision = strsplit(parent,'_');
parentdecision = parentdecision(1);

%Add the parent's decision to the list of previous decisions
previousdecisions = [ListNodes.(parent).previousdecisions, {parentdecision}]; 
previousdecisions = [previousdecisions{:}];

%Determine the possible decisions & the number of times each target can
%still be visisted
[possibledecisions, visitsleft] = DeterminePossDecisions(Inputs, ListNodes, parent, previousdecisions, targetname, attributes(1));

%Create structure of the new node
newNode = struct('node_ID',           node_ID,...                   % The ID of the node
                 'parent',            parent, ...                   % The parent of the node
                 'children',          [],...                        % Matrix that holds the nodes' connections to each other
                 'radius',            [Inputs.StartingRadius],...   % The radius of each connection
                 'length',            [veinlength],...                        % The length of each connection
                 'flux',              [],...                        % Matrix containing each connection's flux
                 'attributes',        [Attributes], ... % Attributes that describe this node (such as orbital elements & ToF .)
                 'previousdecisions', {previousdecisions},...       %List of previous decisions made                    
                 'possibledecisions', {possibledecisions}, ...      %List containing the decisions that can still be made
                 'VisitsLeft',        {visitsleft} ...              % Vector containing the number of times each target cna still be visisted
                ); 

 
%Calculate the flux and add it to the structure
newNode.flux = [CalculateFlux(Inputs,newNode)];
                 
end


