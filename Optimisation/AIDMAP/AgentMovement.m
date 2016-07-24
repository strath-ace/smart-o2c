function [Solutions, ListNodes, Agents, agentdeathflag, funccalls] = AgentMovement(Inputs, Solutions, ListNodes, Agents, agent, funccalls)
%% AgentMovement: This function is used to move agents through the graph
% 
%% Inputs:
% * Inputs      : Structure containing the PhysarumSolver inputs
% * Solutions   : Structure containing the solutions found so far
% * ListNodes   : Structure containing the graph
% * Agents      : Structure containing the Agents
% * agent       : Cell with the current agents' name
% * funccalls   : The number of cost function calls performed so far
% 
%% Outputs:
% * Solutions      : Updated structure containing the solutions found so far
% * ListNodes      : Structure containing the graph
% * Agents         : The structure containing the agents with their new positions
% * agentdeathflag : Flag that shows whether the agent has 'died'
% * funccalls      : The number of cost function calls performed so far
% 
%% Author: Aram Vroom (2016)
% Email: aram.vroom@strath.ac.uk

% Initialise the agent death flag
agentdeathflag = 0;

% For ease of reading, define the current agent and the current node as
% variables
currentagent = char(agent);
currentnode = char(Agents.(currentagent).currentNode);

% Display the agent's movement in the command window
disp(strcat([datestr(now), ' === Moved to', ' ', currentnode]));

% Check if end conditions reached
continueflag = Inputs.EndConditionsFile(Inputs, Agents, agent);

% Stop agent if end conditions reached or no more possible cities to visit
if ((isempty(ListNodes.(currentnode).possibledecisions)) || (sum(ListNodes.(currentnode).VisitsLeft) == 0) || continueflag == 0)
       
    % Set the agent deathflag to 1
    agentdeathflag = 1;
    
    % Save the solution
    Solutions.Nodes = [Solutions.Nodes; {[Agents.(char(agent)).previousListNodes {Agents.(char(agent)).currentNode}]}];   
    Solutions.Costs = [Solutions.Costs; {[Agents.(currentagent).previouscosts]}];
    return
end

% Generate a random number and retreive the ramification probability
p = rand;
pram = Inputs.RamificationProbability;

% Generate a set of potential nodes to ramificate to
[ListNodes, GeneratedNodes, agentdeathflag, funccalls] = Ramification(Inputs, ListNodes, Agents, currentagent, funccalls);

% Retrieve the node IDs of the generated nodes and Initialise a matrix to
% retrieve the fluxes
generatednodenames = fieldnames(GeneratedNodes);
ramfluxesmod = zeros(1, length(generatednodenames));

% Obtain the fluxes of the generated nodes and include the ramification weight
for i = 1:length(generatednodenames)
    ramfluxesmod(i) = GeneratedNodes.(generatednodenames{i}).flux^(Inputs.RamificationWeight);
end
    
% Check whether the current node has children. 
if (~isempty(ListNodes.(currentnode).children))
    
    % Pre-allocation
    childfluxes = zeros(1, length(ListNodes.(currentnode).children));
    
    % Obtain the fluxes of the children
    for i = 1:length(ListNodes.(currentnode).children)
        childfluxes(i) = ListNodes.(char(ListNodes.(currentnode).children(i))).flux;
    end
    
    % Calculate the probability of going to each node. This vector has length 
    % number of childs + number of generated nodes; it contains the probability 
    % for each node. The smaller the cost (higher flux) the higher the probability
    problist = [pram.*ramfluxesmod (1-pram).*childfluxes]./(sum(pram.*ramfluxesmod)+sum((1-pram).*childfluxes));
    
else
    
    % If there are no children, only use the fluxes of the nodes generated
    % through the ramification process
    problist = ramfluxesmod./sum(ramfluxesmod);
end

% Sanity check #2
if isempty(problist)
    
    % If there are no probabilities to choose from (no children and no 
    % generated nodes), set the death flag to 1
    agentdeathflag = 1;
    
    % Save the solution
    Solutions.Nodes = [Solutions.Nodes; {[Agents.(char(agent)).previousListNodes {Agents.((char(agent))).currentNode}]}];   
    Solutions.Costs = [Solutions.Costs; {[Agents.(currentagent).previouscosts]}];
    return
end

% Choose a node and find its index
chosenindex = find(rand<cumsum(problist), 1, 'first');

% If the index is smaller than the number of generated nodes, it falls
% within that list
if chosenindex <=length(generatednodenames)
    
    % Retrieve the ID of the chosen node
    chosennode = generatednodenames{chosenindex};
    
    % Add chosen node to the ListNodes structure
    chosennodestruct = GeneratedNodes.(chosennode);
    ListNodes = AddNewNode(ListNodes, chosennodestruct);
    
% If the index falls outside of the list of generated nodes, the agent
% moves to one of the children    
else
    
    % Retrieve the children of the node
    nodechildren = cell(ListNodes.(currentnode).children);
    
    % Find the child that was chosen. The -length(generatednodenames) is
    % used to find the correct index of the child in the children list of the node 
    chosennode = nodechildren{chosenindex-length(generatednodenames)};

end
  
% Move the agent to the chosen node
Agents.(currentagent).previousListNodes = [Agents.(currentagent).previousListNodes {currentnode}];
Agents.(currentagent).currentNode = chosennode;
Agents.(currentagent).previouscosts = [Agents.(currentagent).previouscosts ListNodes.(chosennode).length];

end
   