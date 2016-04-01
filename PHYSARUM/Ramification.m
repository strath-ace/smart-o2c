function [Nodes,Agents] = Ramification(Inputs,Nodes,Agents,agent)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if isempty(Nodes.ListNodes.(char((Agents.(agent).currentNode))).possibledecisions)
    return
end

%Set boundaries
mincharboundary  = 0;
maxcharboundary = 9;
stepsize = 1;

%Create a list of all possible characteristics & decisions
possiblecharacteristics = mincharboundary:stepsize:maxcharboundary;
possibledecisions = Nodes.ListNodes.(char(Agents.(agent).currentNode)).possibledecisions;

%Create list of nodes that can be selected
for i = 1:(maxcharboundary-mincharboundary+1)
    for j = 1:length(possibledecisions)
    possdeccharvec(i,j) = strcat(possibledecisions(j),'_',num2str(possiblecharacteristics(i)));
    end
end
possdeccharvec = reshape(possdeccharvec,[1,numel(possdeccharvec)]);

%Retrieve list of currently existing nodes
existingnodes = fields(Nodes.ListNodes);

%Initialize empty structures to save the generated nodes in
generatednodes = struct([]);
nameslist = [];
problist = [];

%Disable the "Concatenate empty structure" warning
warning('off','MATLAB:catenate:DimensionMismatch');

%Start loop to generate the desired number of nodes
while (length(fields(generatednodes)) < Inputs.RamificationAmount)
    
    %If no more decisions are possible, exit while loop
    if isempty(possdeccharvec)
        break
    end
       
    %Choose a random decision (node_ID)
    randdecision = randi([1 length(possdeccharvec)]);
    newnode_ID = char(possdeccharvec(randdecision));
    
    %Remove chosen decision from list of possible decisions
    possdeccharvec(strmatch(newnode_ID,possdeccharvec)) = [];

    %Confirm that node doesn't already exist
    if (isempty(strmatch(newnode_ID,existingnodes,'exact')) && isempty(strmatch(newnode_ID,fields(generatednodes),'exact')))
        
        %Split the newnode_ID into the chosen target & characteristic
        temp = strsplit(newnode_ID,'_');
        explosiondecision = char(temp(1)); randchar = str2num(char(temp(2)));
        
        %Generate the new node
        [newNode] = CreateNode(Inputs,Nodes,explosiondecision,randchar,char(Agents.(agent).currentNode));
        newNode.(char(newnode_ID)).lengths = CostFunction(Nodes.ListNodes.(char(Agents.(agent).currentNode)),newNode.(char(newnode_ID)));
        
        %Add generated node to the structure created earlier
        tmp = [fieldnames(generatednodes), struct2cell(generatednodes); fieldnames(newNode), struct2cell(newNode)].';
        generatednodes = struct(tmp{:});
        
        %Add generated node name & cost to matrices for ease of access
        nameslist = [nameslist {newnode_ID}];
        
        %Prevent 1 / 0
        if (newNode.(char(newnode_ID)).lengths == 0)
            newNode.(char(newnode_ID)).lengths = 1e-15;
        end
        
        %Create a vector with the cost to each generated node
        problist = [problist newNode.(char(newnode_ID)).lengths];
    end
end

%If no nodes are found, exit the function
if isempty(problist)
    return
end

%Use node name & cost matrices to make a probabilistic choice based on the
%cost
problist = 1./(problist);
problist = problist./(sum(problist));
chosennode = nameslist(find(rand<cumsum(problist),1,'first'));

%Add chosen node to the Nodes structure
chosennodestruct = struct(char(chosennode),struct(generatednodes.(char(chosennode))));
Nodes = AddNode(Nodes,chosennodestruct);
    
%Move agent to the new node
Agents.(agent).previousNodes = [Agents.(agent).previousNodes Agents.(agent).currentNode];
Agents.(agent).currentNode = chosennode;
    

end



