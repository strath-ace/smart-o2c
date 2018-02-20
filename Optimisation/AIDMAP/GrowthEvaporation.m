function [ListNodes, BestSolution] = GrowthEvaporation(Inputs, ListNodes, Solutions, BestSolution)
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
<<<<<<< HEAD
%-----------Copyright (C) 2018 University of Strathclyde and Authors-----------
=======
%-----------Copyright (C) 2016 University of Strathclyde-------------
>>>>>>> 5b7361d93c9119cf1d2e9e6c885bed93f924d71b
%
%
%
%% GrowthEvaporation: This function finds the best chain of veins used by the agents and simulates the growth factor and evaporation
% 
%% Inputs:
% * Inputs         : Structure containing the PhysarumSolver inputs
% * ListNodes      : Structure containing the graph
% * Solutions      : Structure containing the solutions so far
% * BestSolution   : Structure containing the best solution found
%                  * BestSolution.BestChain = cell array containing
%                     the node IDs of the best chain
%                   * BestSolution.BestCost = the total cost
%                     corresponding to the best chain        
% 
%% Outputs: 
% * ListNodes       : Structure containing the updated graph
% * BestSolution    : Structure containing the best solution found
%                    * BestSolution.BestChain = cell array containing
%                      the node IDs of the best chain
%                    * BestSolution.BestCost = the total cost
%                      corresponding to the best chain      
% 
%% Author(s): Aram Vroom (2016)
% Email:  aram.vroom@strath.ac.uk


%% Growth Factor

% Loop over the solutions and find the number of nodes and costs
numberofnodes = zeros(1, length(Solutions.Nodes));
costs = zeros(1, length(Solutions.Nodes));
for i = 1:length(Solutions.Nodes)
    numberofnodes(i) = length(Solutions.Nodes{i});
    costs(i) = sum(Solutions.Costs{i});
end

% Retrieve the best chain of nodes
[bestchainindex, bestcost] = Inputs.BestChainFile(numberofnodes, costs);

% In case multiple best solutions are found, only take the first one
bestchainindex = bestchainindex(1);
bestcost = bestcost(1);

% Save best solution found
BestSolution.BestChain = Solutions.Nodes{bestchainindex};
BestSolution.BestCost = bestcost;

% Print the best solution so far
disp('-------------------------------------------------------')
disp(['The best chain so far has ', num2str(length(BestSolution.BestChain)), ' nodes and a cost of ', num2str(bestcost(1))]);
disp('-------------------------------------------------------')

% For ease of reading, rename the best chain into another variable    
bestchain = BestSolution.BestChain;

% Loop over the nodes in the best chain
for i = 1:length(bestchain)

    % For ease of reading, define the node currently being evaluated & its parent
    % as a separate variable
    evaluatednode = char(bestchain(i));
    parent = ListNodes.(evaluatednode).parent;

    % Check if the node has a parent - Ignores root
    if ~isempty(parent)

        % Dilate the respective link in the parent's node structure
        ListNodes.(evaluatednode).radius = ListNodes.(evaluatednode).radius + Inputs.GrowthFactor*ListNodes.(evaluatednode).radius;

        % Check if the link's radius is not too large or small. Correct if
        % so, set to min/max radius
        ListNodes.(evaluatednode).radius(ListNodes.(evaluatednode).radius./Inputs.StartingRadius > Inputs.MaximumRadiusRatio) = Inputs.MaximumRadiusRatio*Inputs.StartingRadius;
        ListNodes.(evaluatednode).radius(ListNodes.(evaluatednode).radius./Inputs.StartingRadius < Inputs.MinimumRadiusRatio) = Inputs.MinimumRadiusRatio*Inputs.StartingRadius;
    end

end



%% Evaporation
allnodes = fieldnames(ListNodes);

for i = 1:length(allnodes)
    % For ease of reading, define the node currently being evaluated & its parent
    % as a separate variable
    evaluatednode = char(allnodes(i));
    parent = ListNodes.(evaluatednode).parent;
    if ~isempty(parent)
        % Simulate evaporation in this link    
        ListNodes.(evaluatednode).radius = ListNodes.(evaluatednode).radius - Inputs.EvaporationCoefficient*ListNodes.(evaluatednode).radius;
        
        % Check if the link's radius is not too large or small. Correct if
        % so, set to max/min radius
        ListNodes.(evaluatednode).radius(ListNodes.(evaluatednode).radius./Inputs.StartingRadius > Inputs.MaximumRadiusRatio) = Inputs.MaximumRadiusRatio*Inputs.StartingRadius;
        ListNodes.(evaluatednode).radius(ListNodes.(evaluatednode).radius./Inputs.StartingRadius < Inputs.MinimumRadiusRatio) = Inputs.MinimumRadiusRatio*Inputs.StartingRadius;
        
        % Update flux
        ListNodes.(evaluatednode).flux = CalculateFlux(Inputs, ListNodes.(evaluatednode));
    end
end

end
