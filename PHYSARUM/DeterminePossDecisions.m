function [possibledecisions,visitsleft] = DeterminePossDecisions(Inputs,Nodes,parent,previousdecisions,decisionname)
%This function creates the structure for a new node.
%
% Inputs:
% * Inputs         : Structure containing the PhysarumSolver inputs
% * Nodes          : Structure that contains the Graph's information
% * parent         : a string containing the node_ID of the parent node
% * previousdecisions : a list of the previous decisions made
% * decisionname   : string with the chosen target
%
% Outputs: 
% * possibledecisions : a list of the decisions that can still be made by
%                       the node
% * visitsleft     : updated vector with the number of times each target
%                    can still be visisted
%
% Author: Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk


%Create a concatenated list of the previous decisions, including the newly generated node 
previousdecisions = strcat(previousdecisions{:}, decisionname);

%Get the decision's (target's) maximum consecutive resonant orbits
decisionindex = strmatch(decisionname,Inputs.PossibleDecisions);

%In case the user accidentally has same target multiple times in the list of
%poss decisions, only keep the index of the first 
decisionindex = decisionindex(1);

maxconsecutive = Inputs.MaxConsecutiveRes(decisionindex);

%Copy the number of visits left of the parent
visitsleft = Nodes.(parent).VisitsLeft;

%Find the number of times a target can still be visited
visitsleft(decisionindex) = Nodes.(parent).VisitsLeft(decisionindex)-1;

%Retrieve all the possible decisions
alldecisions = Inputs.PossibleDecisions;

%Initially, set the decisions the node can make equal to be equal to all the decisions
possibledecisions = alldecisions;

%Remove decisions with visitsleft = 0
possibledecisions(visitsleft == 0) = [];

%Check whether resonance orbit is allowed by first generating the string that
%should be checked for in the previous decisions
checkforcontargets=repmat(decisionname,1,maxconsecutive+1);

%See if max. resonant orbit has been performed as last trajectory so far 
resperformed = max(strfind(previousdecisions,checkforcontargets))==(length(previousdecisions)-length(checkforcontargets)+1);

%If so, exclude the respective target from the list of possible decisions
if resperformed 
    decisionindex = strmatch(decisionname,possibledecisions);
    possibledecisions(decisionindex) = [];
end



end

