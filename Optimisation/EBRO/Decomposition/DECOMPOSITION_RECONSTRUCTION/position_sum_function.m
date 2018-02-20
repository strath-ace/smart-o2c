% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [PS] = position_sum_function(num_combination_max, num_functions_tot, num_combination, decomposition_ii, num_sample)


A =1;

for index_num_Belief = num_functions_tot:-1:2
    if isempty(decomposition_ii{num_sample}.step_one{index_num_Belief}) == 0
        
        A = A*length(decomposition_ii{num_sample}.step_one{index_num_Belief});
        PS(index_num_Belief) = ceil(num_combination*A/num_combination_max);
        num_combination = num_combination-num_combination_max/A*(PS(index_num_Belief)-1);
       
    end
end
PS(1) = num_combination;



n_sub = 0;
for i=1:num_functions_tot
    if isempty(decomposition_ii{num_sample}.step_one{i}) == 0
            n_sub = n_sub + 1;
    end
end

PS = PS(PS~=0);
PS = PS(end-n_sub+1 : end);
end