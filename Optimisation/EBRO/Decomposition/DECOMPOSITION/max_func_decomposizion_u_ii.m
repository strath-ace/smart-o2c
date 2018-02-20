% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [decomposition_ii] = max_func_decomposizion_u_ii(i, ii, in, position_FE, decomposition_ii, algo_inner, problem_inner, num_Belief, j, Sample)

global num_maximization_decomposition
num_maximization_decomposition = num_maximization_decomposition +1;

%% MAX
problem_inner.objfun = @mask_objfun_max_decomposition;

[ u_max_to_opt, fmax_to_opt , ~ , output_aux_to_opt ] = algo_inner.optimise(problem_inner,algo_inner.par);


decomposition_ii{num_Belief}.sample{1,j}.n_sample = j;
decomposition_ii{num_Belief}.sample{1,j}.componente{1,i}.function = i; 
decomposition_ii{num_Belief}.sample{1,j}.componente{1,i}.FocalElement{1,ii}.upper_u = u_max_to_opt;           %[decomposition{1,i-3}.FocalElement{1,ii}.ub; problem_inner.par_objfun.umax];
decomposition_ii{num_Belief}.sample{1,j}.componente{1,i}.FocalElement{1,ii}.upper_f = -fmax_to_opt;

bpa = 1;
for k = 1:problem_inner.dim
    bpa = bpa*problem_inner.bpa{1,1}{k,1}(position_FE(k)); 
end


decomposition_ii{num_Belief}.sample{1,j}.componente{1,i}.FocalElement{1,ii}.bpa = bpa; 

decomposition_ii{num_Belief}.sample{1,j}.componente{1,i}.FocalElement{1,ii}.n_FE = ii;
decomposition_ii{num_Belief}.sample{1,j}.componente{1,i}.FocalElement{1,ii}.lb = problem_inner.lb;
decomposition_ii{num_Belief}.sample{1,j}.componente{1,i}.FocalElement{1,ii}.ub = problem_inner.ub;
        
%% MIN
problem_inner.objfun = @mask_objfun_min_decomposition;

[ u_min_to_opt, fmin_to_opt , ~ , ~ ] = algo_inner.optimise(problem_inner,algo_inner.par);

decomposition_ii{num_Belief}.sample{1,j}.componente{1,i}.FocalElement{1,ii}.downer_u = u_min_to_opt;           
decomposition_ii{num_Belief}.sample{1,j}.componente{1,i}.FocalElement{1,ii}.downer_f = fmin_to_opt;

end