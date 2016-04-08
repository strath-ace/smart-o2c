function [Nodes,Agents,agentdeathflag] = Ramification(Inputs,Nodes,Agents,agent)
% This function handles the ramification to new nodes. It does so by
% generating a preset number of random nodes and making a probabilistic
% selection based on the cost function.
%
% Inputs:
% * Nodes       : Structure containing the graph
% * Inputs      : Structure containing the PhysarumSolver inputs
% * Agents      : The structure containing the agents
% * agent       : A string containing the name of the current agent
%
% Outputs: 
% * Nodes       : Nodes structure where the radii have been updated with
%                 the dilation and evaporation
% * Agents      : The updated structure containing the agents
%
% Author: Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk

%If there are no possible decisions, exit function

%For easy of reason, save the agent's current node in a variable
currentNode = char(Agents.(agent).currentNode);
agentdeathflag = 0;

if isempty(Nodes.(currentNode).possibledecisions)   
    agentdeathflag = 1;
    return
end

%Find all the nodes that can be chosen:

%To do so, get first all the possible nodes that can be made over the
%entire graph. Split these nodes into their target & characteristic
possnodes = Inputs.PossibleNodes;
temp = regexp(possnodes, '_', 'split');
[temp]=cat(1,temp{:});

%Next, retrieve the decisions possible in this node
possdecisions = Nodes.(currentNode).possibledecisions;

%Remove all the nodes that do not have as decision one of possible decisions
%in this node
possnodes(ismember(temp(:,1),possdecisions)==0) = [];

%Retrieve list of currently existing nodes
existingnodes = fields(Nodes);

%Initialize structures to save the generated nodes in. The generatednodes
%structure has a temporary field to circumvent issues with adding fields to
%empty structures
generatednodes = struct('temp',0);
nameslist = [];
costvec = [];

%Disable the "Concatenate empty structure" warning
warning('off','MATLAB:catenate:DimensionMismatch');

%Start loop to generate the desired number of nodes
while (length(fields(generatednodes)) < Inputs.RamificationAmount)
    
    %If no more decisions are possible, exit while loop and set
    %agentdeathflag to 1
    if isempty(possnodes)
        disp(strcat(agent,' died'))
        agentdeathflag = 1;
        break
    end
       
    %Choose a random decision (node_ID)
    randdecision = randi([1 length(possnodes)]);
    newnode_ID = char(possnodes(randdecision));
    
    %Remove chosen decision from list of possible decisions
    possnodes(strmatch(newnode_ID,possnodes)) = [];

    %Confirm that node doesn't already exist
    if (isempty(strmatch(newnode_ID,existingnodes,'exact')) && isempty(strmatch(newnode_ID,fields(generatednodes),'exact')))
        
        %Split the newnode_ID into the chosen target & characteristic
        temp = strsplit(newnode_ID,'_');
        explosiondecision = char(temp(1)); randchar = str2num(char(temp(2)));
        
        %Generate the new node & save its cost in a vector
        [newNode] = CreateNode(Inputs,Nodes,explosiondecision,randchar,currentNode);
        costvec = [costvec CostFunction(Nodes.(currentNode),newNode)];
        
        %Add generated node to the structure created earlier.
        generatednodes.(newNode.node_ID) = newNode;
        
        %Add generated node name & cost to matrices for ease of access
        nameslist = [nameslist {newnode_ID}];
        
        %Prevent 1 / 0
        if (costvec(end) == 0)
            costvec(end) = 1e-15;
        end      
    end
end

%Remove temporary field within the generatednodes stucture
generatednodes = rmfield(generatednodes,'temp');

%If no nodes are found, exit the function
if isempty(costvec)
    return
end

%Use node name & cost matrices to make a probabilistic choice based on the
%cost. The lower the cost, the higher the chance of the node being selected
problist = 1./(costvec.^Inputs.RamificationWeight);
problist = problist./sum(problist);

chosenindex = find(rand<cumsum(problist),1,'first');
chosennode = char(nameslist(chosenindex));

%Add chosen node to the Nodes structure
chosennodestruct = generatednodes.(chosennode);
Nodes = AddNode(Inputs,Nodes,chosennodestruct,costvec(chosenindex));

%Move agent to the new node
Agents.(agent).previousNodes = [Agents.(agent).previousNodes {currentNode}];
Agents.(agent).currentNode = chosennode;
Agents.(agent).previouscosts = [Agents.(agent).previouscosts costvec(chosenindex)];


end


