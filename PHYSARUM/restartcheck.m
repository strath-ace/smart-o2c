function [restartflag] = restartcheck(ListNodes,currentnode)
% This function handles the restarting of the Physarum algorithm
%
% Inputs:
% * ListNodes       : Structure containing the graph
% * currentnode : The current node the agent is at
%
% Outputs: 
% * restartflag : A flag that is set to 1 if the algorithm is to be
%                 restarted
%
% Author: Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk

%Initialize the restart flag
restartflag = 0;

%Check if no more decisions are available
if isempty(ListNodes.(char(currentnode)).possibledecisions)
    
    %Set the restart flag to 1
    restartflag = 1;
    return
end


end