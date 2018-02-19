% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [surrogate] = surrogate_update(surrogate)

for i_FE = 1:size(surrogate.model, 2)
    if surrogate.changed_FE(i_FE) > 0
        % filter the DOE
        [~, idx] = unique(round(1e8*surrogate.x_doe{i_FE}),'rows');
        surrogate.x_doe{i_FE} = surrogate.x_doe{i_FE}(idx,:);
        surrogate.f_doe{i_FE} = surrogate.f_doe{i_FE}(idx,:);
        % retrain the model
        surrogate.model{i_FE} = surrogate.training(surrogate.x_doe{i_FE},surrogate.f_doe{i_FE},surrogate);       
        surrogate.changed_FE(i_FE) = 0;
    end
end