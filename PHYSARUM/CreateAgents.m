function [Agents] = CreateAgents(NumberOfAgents)
%This function creates the structure containing the Agent decisions and
%other characteristics.
%
% Inputs:
% * NumberOfAgents : the number of agents used in the solver
%
% Outputs: 
% * Agents         : the structure containing the set of agents and their
%                    characteristics
%
% Author: Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk

%Initialize the Agents structure
Agents = struct();          

%Loop over the number of agents to be created
for i = 1:NumberOfAgents
    
    %Find the field name of the specific agent in the Agents structure
    strcurrid = ['agent',num2str(i)];
    
    %Create the structure for each agent
    AgentStruct = struct(strcurrid, struct(               ...
                                      'previousNodes',[], ...
                                      'previouscosts',[], ...
                                      'currentNode',  ['Root'],  ...
                                      'currentconnections',[]...
                                      ));
    
    %Add the structure for this agent underneath the structures for the
    %previously generated agents
    tmp = [fieldnames(Agents), struct2cell(Agents); fieldnames(AgentStruct), struct2cell(AgentStruct)].';
    Agents = struct(tmp{:});
end

