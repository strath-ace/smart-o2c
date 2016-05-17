function [ListNodes] = CreateListNodes(Inputs)
% This function creates the structure describing the graph.
%
% Inputs:
% * Inputs : Structure containing the PhysarumSolver inputs
%
% Outputs: 
% * ListNodes  : Empty structure that will contain the Graph's information
%
% Author: Aram Vroom - 2016
% Email:  aram.vroom@strath.ac.uk


%%%Create UID%%%
%As only underscores can be used in field names, the UID is made up as follows:
%3 underscores define the difference between the target and the
%  attributes
%2 underscores define the difference between two attributes
%1 underscore denotes the location of the decimal point within an
%  attribute

rootattribstr = [];
for i = 1:(length(Inputs.RootAttrib)*2-1)
    if rem(i,2)
        rootattribstr = strcat(rootattribstr,num2str(Inputs.RootAttrib((i+1)/2)));
    else
        rootattribstr = strcat(rootattribstr,'__');
    end
end

rootattribstr = strrep(rootattribstr,'.','_');

rootID = strcat(Inputs.RootName,'___',rootattribstr);

%Create the ListNodes structure
ListNodes = struct(struct(rootID,   ...
                         struct(...
                         'node_ID',           rootID,...
                         'parent',            [],... % The parent of the node
                         'children',          [],... % Matrix that holds the nodes' connections to each other
                         'radius',            [],... % The radius of each connection
                         'length',           [],... % The length of each connection
                         'flux',            [],... % Matrix containing each connection's flux
                         'attributes',      [SetNodeAttributes(Inputs,[],Inputs.RootName,Inputs.RootAttrib)],... % Attributes that describe this node (such as orbital elements & ToF .)
                         'previousdecisions', [],... % List of the previous decisions made
                         'possibledecisions', {Inputs.PossibleDecisions}, ... % Targets that can still be visisted by the node
                         'VisitsLeft',        {Inputs.MaxVisits} ... % Vector containing the number of times each target cna still be visisted
                     )));
           
end 
