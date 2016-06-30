function [generatednodes, agentdeathflag, funccalls] = Ramification2(Inputs, Solutions, ListNodes, Agents, agent, funccalls)
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

%For easy of reading, save the agent's current node in a variable
currentNode = char(Agents.(agent).currentNode);
agentdeathflag = 0;

%If there are no possible decisions, exit function
if (isempty(ListNodes.(currentNode).possibledecisions))
    agentdeathflag = 1;
    return
end

%Find all the nodes that can be chosen:

%To do so, get first all the IDs that can be in the tree
temp = strsplit(currentNode,'___');
possnodes = Inputs.PossibleListNodes;
possids = strcat(temp(end),'___',possnodes);

%Create a vector that tracks the index of the possible nodes and their
%feasibility. If a node is not feasible, its corresponding index in this
%vector will be set to NaN
indextracker = 1:length(possnodes);

%Set the index of the already existing nodes to NaN
indextracker(ismember(possids,fieldnames(ListNodes)))=NaN;

%Split the remaining nodes into their target & characteristic
temp = regexp(possnodes, '__', 'split');
[temp]=cat(1, temp{:});

%Next, retrieve the decisions possible in this node
possdecisions = ListNodes.(currentNode).possibledecisions;

%Set the indices to NaN of the nodes that do not have as decision one of 
%possible decisions in this node
indextracker(ismember(temp(:,1), possdecisions)==0) = NaN;

%Initialize structures to save the generated nodes in. The generatednodes
%structure has a temporary field to circumvent issues with adding fields to
%empty structures
generatednodes = struct('temp',0);
nameslist = cell(1,Inputs.RamificationAmount);

%Initial index for the nameslist variable, the attempt counter and a
%tracker for the number of NaNs in the indextracker vector
i = 1;
attempt = 1;
nantracker = sum(isnan(indextracker));

%Disable the "Concatenate empty structure" warning
warning('off','MATLAB:catenate:DimensionMismatch');

%Start loop to generate the desired number of nodes
while (length(fields(generatednodes)) <= Inputs.RamificationAmount)

    %If all values in indextracker are NaN (meaning no more possible
    %children to choose from) or max number of attempts has been reached, exit while loop
    if ((nantracker == length(indextracker)) || (attempt == Inputs.MaxChildFindAttempts))      
        break
    end
    
    %Choose a node from the list of possible nodes to generate
    [newnode_ID,nodeindex] = ChooseNode(Inputs,currentNode,possnodes,indextracker,nantracker);
    
    %Increase the attempt & NaN counter by 1
    attempt = attempt +1;
    nantracker = nantracker+1;

    %Remove chosen decision from list of possible decisions
    indextracker(nodeindex) = NaN;
    
    %Check if the node is valid based on the UID
    [validflag] = Inputs.NodeIDCheckFile(Inputs,ListNodes,newnode_ID,currentNode,generatednodes);
   
    %Confirm that node doesn't already exist
    if (~validflag)
        continue
    end
    
    funccalls = funccalls+1;
    
    %Generate the new node & save its cost in a vector
    [newNode] = CreateNode(Inputs, ListNodes, newnode_ID, currentNode);
    
    
    if isempty(newNode)
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



