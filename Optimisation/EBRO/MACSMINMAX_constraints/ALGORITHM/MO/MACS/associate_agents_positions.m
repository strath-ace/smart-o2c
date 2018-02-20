% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [id_point,id_agent_out,vals] = associate_agents_positions (g_vals,id_agent_in)

n_iter = min(size(g_vals));
gtmp_star = g_vals;

id_point = zeros(n_iter,1);
id_agent_out = zeros(n_iter,1);
vals = zeros(n_iter,1);

for i = 1:n_iter

    [a,b] = min(gtmp_star(:));
    [I_row, I_col] = ind2sub(size(gtmp_star),b);    % locate best association
    
    id_point(i) = I_row;                            % associate point to problem
    id_agent_out(i) = id_agent_in(I_col);
    vals(i) = a;
    
    gtmp_star(I_row,:) = Inf;                        % trick to remove row and column from matrix without actually resizing it
    gtmp_star(:,I_col) = Inf;
    
end

end