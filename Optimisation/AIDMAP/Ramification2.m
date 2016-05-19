function [generatednodes, agentdeathflag] = Ramification2(Inputs, Solutions, ListNodes, Agents, agent)
% This function handles the ramification to new nodes. It does so by
% generating a preset number of random nodes and making a probabilistic
% selection based on the cost function.
%
% Inputs:
% * ListNodes       : Structure containing the graph
% * Inputs      : Structure containing the PhysarumSolver inputs
% * Agents      : The structure containing the agents
% * agent       : A string containing the name of the current agent
%
% Outputs: 
% * ListNodes       : ListNodes structure where the radii have been updated with
%                 the dilation and evaporation
% * Agents      : The updated structure containing the agents
%
% Author: Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk

%If there are no possible decisions, exit function

%For easy of reading, save the agent's current node in a variable
currentNode = char(Agents.(agent).currentNode);
agentdeathflag = 0;

if (isempty(ListNodes.(currentNode).possibledecisions))
    agentdeathflag = 1;
    return
end

%Find all the nodes that can be chosen:

%To do so, get first all the possible nodes that can be made over the
%entire graph. 
possnodes = Inputs.PossibleListNodes;

%Remove the already existing nodes
possnodes(ismember(possnodes,fieldnames(ListNodes)))=[];

%Split the remaining nodes into their target & characteristic
temp = regexp(possnodes, '___', 'split');
[temp]=cat(1, temp{:});

%Next, retrieve the decisions possible in this node
possdecisions = ListNodes.(currentNode).possibledecisions;

%Remove all the nodes that do not have as decision one of possible decisions
%in this node
possnodes(ismember(temp(:,1), possdecisions)==0) = [];

%Initialize structures to save the generated nodes in. The generatednodes
%structure has a temporary field to circumvent issues with adding fields to
%empty structures
generatednodes = struct('temp',0);
nameslist = cell(1,Inputs.RamificationAmount);

%Initial index for the nameslist variable and the attempt counter
i = 1;
attempt = 1;

%Disable the "Concatenate empty structure" warning
warning('off','MATLAB:catenate:DimensionMismatch');

%Start loop to generate the desired number of nodes
while (length(fields(generatednodes)) <= Inputs.RamificationAmount)
    
    %If no more decisions are possible, exit while loop and set
    %agentdeathflag to 1
    if isempty(possnodes)
        disp(strcat(agent,' died'))
         agentdeathflag = 1;
        break
    end
    
    if (attempt == Inputs.MaxChildFindAttempts)
        break
    end
    
    %Increase the attempt counter by 1
    attempt = attempt +1;
    
    %Choose a node from the list of possible nodes to generate
    [newnode_ID,childID] = ChooseNode(currentNode,possnodes);
    
    %Remove chosen decision from list of possible decisions
    possnodes(strcmp(childID,possnodes)) = [];
    
    %Check if the node is valid based on the UID
    [validflag] = Inputs.NodeIDCheckFile(ListNodes,newnode_ID,currentNode,generatednodes);
   
    %Confirm that node doesn't already exist
    if (~validflag)
        continue
    end

        
    %Generate the new node & save its cost in a vector
    [newNode] = CreateNode(Inputs, ListNodes, newnode_ID, currentNode);
    
    %Sanity Check
    if newNode.length == Inf
        continue
    end
    
    %If the node is valid according to the first check, use the created
    %node structure to further determine the validity
    [newNode] = Inputs.CreatedNodeCheckFile(Inputs, newNode, ListNodes);
    
    %Sanity Check
    if newNode.length == Inf
        continue
    end

    %Add generated node to the structure created earlier.
    generatednodes.(newNode.node_ID) = newNode;
    
    %Add generated node name & cost to matrices for ease of access
    nameslist{i} = newnode_ID;           
    
    %Increase the index for the nameslist struct
    i = i+1;
end

%Remove temporary field within the generatednodes stucture
generatednodes = rmfield(generatednodes, 'temp');

end



