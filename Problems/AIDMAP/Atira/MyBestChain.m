function [bestchainindex, bestcost] = MyBestChain(numberofnodes, nodecosts)
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

