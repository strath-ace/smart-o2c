function [restartflag] = RestartCheck(Inputs, ListNodes, Agents)
% This function checks the restarting of the Physarum algorithm
%
% Inputs:
% * Inputs   : Structure containing the PhysarumSolver inputs
% * ListNodes       : Structure containing the graph
% * Agents      : Structure containing the Agents & their characteristics
%
% Outputs: 
% * restartflag : A flag that is set to 1 if the algorithm is to be
%                 restarted
%
% Author: Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk

%Initialize the restart flag
restartflag = 0;

%Retrieve the agent names
agents = fieldnames(Agents);

%Loop over all the agents
for i = 1:length(agents)
    
    %Change the agent into characters (needed for structures)
    agent = char(agents(i));
    
    %Retrieve the nodes visisted by the agent
    visistednodes = [Agents.(agent).previousListNodes {Agents.(agent).currentNode}];
    
    %Loop over all the visisted nodes
    for j = 1:length(visistednodes)
        
        %Obtain the decisions made for these nodes
        temp = strsplit(char(visistednodes(j)),'_');
        decisions(i,j) = temp(1);
        
        %Retrieve the characteristic that determines whether two nodes are
        %equal if they have the same target (eg ToA)
        determiningcharacteristic(i,j) = ListNodes.(char(visistednodes(j))).characteristics(Inputs.DeterminingCharacteristic);
        
        %Concatenate the decision & the determining characteristic to
        %create a matrix used to check whether the same sequence was chosen
        equalcheck(i,j) = strcat(decisions(i,j), num2str(determiningcharacteristic(i,j)));
    end
end
%Loop over all the agents
for i = 1:length(agents)
    for j = 1:length(agents)
        
        %Only compare different agents
        if i ~= j
            
            %Find the number of equal nodes in the sequences
            numberequal = sum(strcmp(equalcheck(i,:), equalcheck(j,:)));
            
            %If this exceeds the set threshold, set the restartflag to 1
            if numberequal > Inputs.MinCommonNodesThres
                restartflag = 1;
            end
        end
    end
end

%For debug purposes, show that a restart was confirmed to be needed
if restartflag
    disp('Restart Needed')
end
        
end