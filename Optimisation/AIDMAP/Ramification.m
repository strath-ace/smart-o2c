function [ListNodes, GeneratedNodes, agentdeathflag, funccalls] = Ramification(Inputs, ListNodes, Agents, agent, funccalls)
%% Ramification: This function handles the ramification to new nodes. 
% It does so by generating a preset number of random nodes and making a probabilistic
% selection based on the cost function.
% 
%% Inputs:
% * Inputs      : Structure containing the PhysarumSolver inputs
% * ListNodes   : Structure containing the graph
% * Agents      : The structure containing the agents
% * agent       : A string containing the name of the current agent
% * funccalls   : The number of cost function calls performed so far
% 
%% Outputs: 
% * ListNodes           : ListNodes structure where the radii have been updated with
%                         the dilation and evaporation
% * GeneratedNodes      : Structure containing the ramification nodes generated
% * agentdeathflag      : Indicator as to whether the agent couldn't move.
%                         Is 1 if this is the case, 0 otherwise
% 
%% Author: Aram Vroom (2016)
% Email:  aram.vroom@strath.ac.uk

% For easy of reading, save the agent's current node in a variable
currentNode = char(Agents.(agent).currentNode);
agentdeathflag = 0;

% If there are no possible decisions, exit function
if (isempty(ListNodes.(currentNode).possibledecisions))
    agentdeathflag = 1;
    return
end

% Find all the nodes that can be chosen:

% To do so, get first all the IDs that can be in the tree
temp = strsplit(currentNode, '___');
possnodes = Inputs.PossibleListNodes;
possids = strcat(temp(end), '___', possnodes);

if Inputs.LowMem == 1
    indextracker = 1:length(possnodes);
else
    indextracker = ListNodes.(currentNode).ChildValidityTracker;
end

% Set the index of the already existing nodes to NaN
indextracker(ismember(possids, fieldnames(ListNodes)))=NaN;

% Split the remaining nodes into their city & characteristic
temp = regexp(possnodes, '__', 'split');
[temp]=cat(1, temp{:});

% Next, retrieve the decisions possible in this node
possdecisions = ListNodes.(currentNode).possibledecisions;

% Set the indices to NaN of the nodes that do not have as decision one of 
% possible decisions in this node
indextracker(ismember(temp(:, 1), possdecisions)==0) = NaN;

% Initialise structures to save the generated nodes in. The generatednodes
% structure has a temporary field to circumvent issues with adding fields to
% empty structures
GeneratedNodes = struct('temp', 0);
nameslist = cell(1, Inputs.RamificationAmount);

% Initial index for the nameslist variable, the attempt counter and a
% tracker for the number of NaNs in the ListNodes.(currentNode).ChildValidityTracker vector
i = 1;
attempt = 1;
nantracker = sum(isnan(indextracker));

% Disable the "Concatenate empty structure" warning
warning('off', 'MATLAB:catenate:DimensionMismatch');

% Start loop to generate the desired number of nodes
while (length(fields(GeneratedNodes)) <= Inputs.RamificationAmount)

    % If all values in ListNodes.(currentNode).ChildValidityTracker are NaN (meaning no more possible
    % children to choose from), a max number of attempts has been reached or
    % all the possible nodes have been generated, exit while loop. The -1 in
    % the last check is due to the temporary field in the generatednodes
    % structure
    if ((nantracker == length(indextracker)) || (attempt == Inputs.MaxChildFindAttempts) || (length(fieldnames(GeneratedNodes))-1) == (length(indextracker)-nantracker))      
        break
    end
    
    % Choose a node from the list of possible nodes to generate
    [newnode_ID, nodeindex] = ChooseNode(Inputs, currentNode, possnodes, indextracker, nantracker);
        
    % Increase the attempt counter by 1
    attempt = attempt +1;
    
    % If the chosen node is already part of the generatednodes structure, 
    % continue to another attempt
    if any(strcmp(fieldnames(GeneratedNodes), newnode_ID))
        continue
    end   
       
    
    % Check if the node is valid based on the UID
    [validflag] = Inputs.NodeIDCheckFile(Inputs, ListNodes, newnode_ID, currentNode);
   
    % Confirm that node doesn't already exist
    if (~validflag)
        % Remove chosen decision from list of possible decisions & increment
        % nantracker
        indextracker(nodeindex) = NaN;
        nantracker = nantracker+1;
        continue
    end
    
    % Update the number of function calls
    funccalls = funccalls+1;
    
    % Generate the new node & save its cost in a vector
    [newNode] = CreateNode(Inputs, ListNodes, newnode_ID, currentNode);
    
        
    if isempty(newNode)
        % Remove chosen decision from list of possible decisions & increment
        % nantracker
        indextracker(nodeindex) = NaN;
        nantracker = nantracker+1;
        continue
    end

    % Add generated node to the structure created earlier.
    GeneratedNodes.(newNode.node_ID) = newNode;
    
    % Add generated node name & cost to matrices for ease of access
    nameslist{i} = newnode_ID;           
    
    % Increase the index for the nameslist struct
    i = i+1;

end

% Remove temporary field within the generatednodes stucture
GeneratedNodes = rmfield(GeneratedNodes, 'temp');

% Save indextracker to the ChildValidityTracker of the node, so the next
% agent can use this information. The code has been written such that this
% does not change the probability for this other agent to find valid nodes
% (as this would otherwise mean information is shared where the algorithm
% should not); saving this information merely saves on computation time
if Inputs.LowMem == 0
    ListNodes.(currentNode).ChildValidityTracker = indextracker;
end

end



