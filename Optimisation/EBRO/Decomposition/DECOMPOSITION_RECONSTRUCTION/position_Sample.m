% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [PS] = position_Sample(num_sample_tot, num_coupled_vectors, num_sample, Sample, in)
%
% 
% find the vector of the position of an element in a set: 
% for an element N in a set M elements divided in S sub-sets is:
%
% position_Sample = [N_1 N_2 ... N_S]
%
% INPUT: - num_sample_tot        : 
%        - num_coupled_vectors
%        - num_sample
%        - Sample        
%        - in
% OUTPUT: PS: vector-position of the sample


A =1;
for num_Belief = num_coupled_vectors:-1:2
    if in.dim_u_i(in.num_functions + num_Belief)>0
    if in.output == 0 || in.output == 2
        A = A*length(Sample{num_Belief}.FE);
    elseif in.output == 1
        A = A*length(Sample{num_Belief}.FE_Plausibility);
    end
    PS(num_Belief) = ceil(num_sample*A/num_sample_tot);
    num_sample = num_sample-num_sample_tot/A*(PS(num_Belief)-1);
    end
end
if in.dim_u_i(in.num_functions + 1)>0
    PS(1) = num_sample;
else
    PS(1) = 0;
end



end

% [i(1),i(2),i(3),i(4),i(5),i(6),i(7),i(8)] = ind2sub(size(max_mat), num_sample);