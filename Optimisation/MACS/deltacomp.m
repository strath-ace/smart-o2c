function [Delta,sD]=deltacomp(A,x,out)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%
% Computes and sorts differences between agent x and all agents in memory
% A. If out is given, element 'out' is penalized, so to make it last after
% sorting.

% INPUT
%       A       :   memory
%       x       :   agent wrt compute differences
%       out     :   OPTIONAL agent to be excluded from memory
%
% OUTPUT
%       Delta   :   vector of distances, sorted
%       sD      :   sorting index, wrt initial indexing
%
% Revised, cleaned and optimized by Lorenzo A. Ricciardi 2015


nA = length(A(:,1));                                                        % number of agents in memory
Delta=zeros(1,nA)+1e137;                                                    % init of vector Delta... 1e137 is a very high number, that eventually will be associated to out-th element if it's given
lx=length(x(1,:));                                                          % number of parameters of agents

Delta(1:end~=out) = (sum((repmat(x,nA-(out~=0),1)-A(1:end~=out,1:lx)).^2,2)).^0.5;  % fast computation of Euclidean distances on all elements different from 'out'.

[Delta,sD]=sort(Delta);                                                     % Distances are sorted, and sD are sorting index wrt initial indexing

return
