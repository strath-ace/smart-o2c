% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [y, mse] = kriging_predictor(x, model)

%wrapper for the predictor function, only for the formatting of outputs
%etc.

    % if iscell(model)
    %     n_obj = length(model);
    %     y=zeros(1,n_obj);
    %     mse = zeros(1,n_obj);
    %     for obj = 1:n_obj
    %         if size(x,1) == 1
    %             [y_aux, ~, mse_aux] = predictor(x,model{obj});
    %             y_aux = y_aux';
    %         else
    %             [y_aux, mse_aux] = predictor(x, model{obj});
    %         end
    %         y(obj)=y_aux;
    %         mse(obj)=mse_aux;
    %     end
    % else
        if size(x,1) == 1
            [y, ~, mse] = predictor(x,model);
            y = y';
        else
            [y, mse] = predictor(x, model);
        end
    % end


return