function [bestchainindex, bestcost] = MyBestChainMainBelt(numberofnodes, nodecosts)
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%-----------Copyright (C) 2016 University of Strathclyde-------------
%
%
%
%% MyBestChainMainBelt: This function finda the best chain to be used in the growth factor functionality
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

% Find the index of the solution with the moster asteroids and the minimum
% cost
bestchainindex = maxasteroidnumindex(minimumcostindex);

% Retreive the best cost
bestcost = maxasteroidcosts(minimumcostindex);

end

