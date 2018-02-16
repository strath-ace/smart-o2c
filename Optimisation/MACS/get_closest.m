function selection = get_closest(targets,candidates)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%
% returns indexes of closest (target/candidate) pairs

if size(targets,2)~=size(candidates,2)
    
    error('Targets and candidates have different dimensions');
    
end

ntarg = size(targets,1);
ncand = size(candidates,1);
ncoord = size(targets,2);

availtarg = ntarg;
availcand = ncand;
availcomb = min(availtarg,availcand);
  
selection = zeros(availcomb,2);

% FAST computation of pairwise distances
xx = permute(targets,[1,3,2]);                                         
xx = repmat(xx,1,size(candidates,1));
yy = permute(candidates,[3,1,2]);
yy = repmat(yy,size(targets,1),1);
M = sum((xx-yy).^2,3);

MM = M;
tidx = 1:ntarg;
cidx = 1:ncand;

for i = 1:availcomb

    % find best within available combinations
    
    [~,idx]=min(MM(:));
    [row,col]=ind2sub(size(MM),idx);

    % associate row/col with real indexes of the targets/candidates
    tchosenindex = tidx(row);
    cchosenindex = cidx(col);
    
    % add choiche to selection
    selection(i,:) = [tchosenindex,cchosenindex];

    % remove used row/col (target/candidate) pairs
    tidx(row) = [];
    cidx(col) = [];
    MM(row,:) = [];
    MM(:,col) = [];
    
end

end