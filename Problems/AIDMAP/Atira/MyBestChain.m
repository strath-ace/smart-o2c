function [bestchainindex, bestcost] = MyBestChain(numberofnodes, nodecosts)
%% MyBestChain: This function finds the best chain to be used in the growth factor functionality
% 
%% Inputs:
% * numberofnodes : Vector with the amount of nodes for each solution so far
% * nodecosts     : Vector with the costs for each solution
% 
%% Outputs: 
% * bestchainindex  : The index of the best chain within the numberofnodes
%                     and nodecosts vectors. Note: AIDMAP does not support 
%                     multiple best chains in a single generation [integer]
% * bestcost        : The cost of the best solution [real number]
% 
%% Author(s): Aram Vroom (2016)
% Email:  aram.vroom@strath.ac.uk

% Find the index of the solution(s) with the most asteroids
maxasteroidnumindex = find(numberofnodes==max(numberofnodes));

% Find the corresponding cost(s)
maxasteroidcosts = nodecosts(maxasteroidnumindex);

% Obtain the index of the minimum cost
minimumcostindex = find(maxasteroidcosts==min(maxasteroidcosts));

% Find the index of the solution with the most asteroids and the minimum
% cost
bestchainindex = maxasteroidnumindex(minimumcostindex);

% Retreive the best cost
bestcost = maxasteroidcosts(minimumcostindex);

end

