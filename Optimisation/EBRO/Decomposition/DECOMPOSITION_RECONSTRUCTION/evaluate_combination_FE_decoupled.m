% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [decomposition_end] = NEW_evaluate_combination_FE_decoupled(decomposition_end, num_sample, problem_decomposition, minmax, minmin, u_max_tot, Sample, PS, u_max_tot_Plausibility, n_obj, n_point, in, Partial_curve)


% FIND number of combinations of Focal Elements of the uncoupled variables
% for each samples of the coupled parameters.
num_combinations_max = 1;
num_functions_tot = length(decomposition_end{num_sample}.step_one);

num_functions_tot_new = 0; % number of non empty subsystems 
for num_functions = 1:num_functions_tot
    if isempty(decomposition_end{num_sample}.step_one{num_functions})==0
        num_combinations_max = num_combinations_max*length(decomposition_end{num_sample}.step_one{num_functions});
        
        num_functions_tot_new = num_functions_tot_new + 1;
        sub_not_empty(num_functions_tot_new) = num_functions;
    end
end
%num_functions_tot = num_functions_tot_new;



decomposition_end{num_sample}.final_FE     = [];
decomposition_end{num_sample}.FE_add_table = [];

%--------------------------------------------------------------------------
% BELIEF
%--------------------------------------------------------------------------
if in.output == 0 || in.output == 2
    
    d_belief = minmax.d;
    f_cost = problem_decomposition.objfun{n_obj}(d_belief, u_max_tot, problem_decomposition);
    
    
    
    % EVALUATE bpa-scale: for each sample in the coupled domain (U_c) do
    % evaluate a scaled curve. The scale factor is the product of the
    % Belief value of each sample in the partial curves
    
    bpa_scale = 1;
    Sample_u_lb = [];
    Sample_u_ub = [];
    
    for j = 1:length(PS)
        if in.dim_u_i(in.num_functions + j)>0
            if PS(j) == 1
                bpa_scale = bpa_scale*Sample{j}.Belief_f(PS(j));
                
            else
                bpa_scale = bpa_scale*(Sample{j}.Belief_f(PS(j)) - Sample{j}.Belief_f(PS(j)-1));
                
            end
            Sample_u_lb = [Sample_u_lb  Sample{j}.FE{PS(j)}.lb];
            Sample_u_ub = [Sample_u_ub  Sample{j}.FE{PS(j)}.ub];
        end
    end
    
    
    % RECONSTRUCTION of all the focal elements of the final scaled curves
    
    % initialise for add-FE
    %     bpa_add = 0;
    F_max = 0;
    F_ij_add2Delta_ij                = [];
    bpa_add2Delta_ij                 = [];
    position_uncoupled__add2Delta_ij = [];
    position_coupled_add2Delta_ij    = [];
    N_add = 1;
    
    % combination of the g_i
    for num_combination = 1:num_combinations_max
        
        [P_FE] = position_sum_function(num_combinations_max, num_functions_tot, num_combination, decomposition_end, num_sample);
        
        bpa =1;
        F = 0;
        upper_u_tot_decoupled = [];
        lower_bound_tot_decoupled = [];
        upper_bound_tot_decoupled = [];
        position_u_uncoupled = [];
        
        for i = 1:length(P_FE)
            bpa = bpa*decomposition_end{num_sample}.step_one{sub_not_empty(i)}{P_FE(i)}.bpa;
            F = F + decomposition_end{num_sample}.step_one{sub_not_empty(i)}{P_FE(i)}.upper_f;
            upper_u_tot_decoupled = [upper_u_tot_decoupled decomposition_end{num_sample}.step_one{sub_not_empty(i)}{P_FE(i)}.upper_u];
            lower_bound_tot_decoupled = [lower_bound_tot_decoupled decomposition_end{num_sample}.step_one{sub_not_empty(i)}{P_FE(i)}.lb];
            upper_bound_tot_decoupled = [upper_bound_tot_decoupled decomposition_end{num_sample}.step_one{sub_not_empty(i)}{P_FE(i)}.ub];
            position_u_uncoupled = [position_u_uncoupled  decomposition_end{num_sample}.step_one{sub_not_empty(i)}{P_FE(i)}.position_FE];
        end
        lower_bound = [lower_bound_tot_decoupled  Sample_u_lb]; %Sample{num_sample}.FE{1, 1}.lb];
        upper_bound = [upper_bound_tot_decoupled  Sample_u_ub]; %Sample{num_sample}.FE{1, 1}.ub];
        
        N=0;
        for i=1:num_functions_tot
            if in.dim_u_i(i)>0
                N = N+1;
            end
        end
        F_recombination_N = F - (N-1)*f_cost;
        decomposition_end{num_sample}.step_two{num_combination}.F = F_recombination_N;
        
        F_belief(num_combination) = decomposition_end{num_sample}.step_two{num_combination}.F;
        
        if F_belief(num_combination) > F_max
            position_F_max = num_combination;
            F_max = decomposition_end{num_sample}.step_two{num_combination}.F;
        end
        
        P_partial = [];
        %         for ii=1:length(decomposition_end{num_sample}.step_one)
        %               P_partial = [P_partial decomposition_end{num_sample}.step_one{ii}{1, 1}.n_FE  Sample{num_sample}.FE{num_sample}];
        for j = 1:length(PS)
            if in.dim_u_i(in.num_functions + j)>0
                P_partial = [P_partial Sample{j}.position(PS(j))];% [P_partial Partial_curve{j}.position_u_coupled{PS(j)}];
            end
        end
        %         end
        
        
        
        
        %% change (?) the 2 "if" condition
        coupled_vectors = in.dim_u_i(in.num_functions +1 : end);
        position_non_zero_c_v = find(coupled_vectors>0);
        num_coupled_vectors = length(coupled_vectors);%(coupled_vectors>0));
        num_comb_add = 1;
        
        for jj = position_non_zero_c_v
            num_comb_add = num_comb_add*length(Partial_curve{jj}.Belief_FE_function);
        end
        %%
        position = [];
        for ii=1:length(Sample)
            if isempty(Sample{ii}) == 0
                position = [position Sample{ii}.position(PS(ii))];
            end
        end
        
        position_next = position;
        for ii=1:length(Sample)
            if isempty(Sample{ii}) == 0
                if Sample{ii}.position(PS(ii))~=Sample{ii}.position(end)
                    position_next = [position_next Sample{ii}.position(PS(ii)+1)];
                end
            end
        end
        
        mm_position = 1;
        nn = 1;
        for mm=1:length(Partial_curve)
            if isempty(Partial_curve{mm}) == 0
                structure{mm}.FE = Partial_curve{mm}.Belief_FE_function;
                mm_position = mm_position +1;
                
                nn = nn + 1;
            end
        end
        
        

        coupled = [];
        for kkk = 1: length(PS)
            if in.dim_u_i(in.num_functions + kkk)>0
%                 coupled = [coupled  Partial_curve{kkk}.u_coupled(PS(kkk),:)];
                    coupled = [coupled  Sample{1, kkk}.FE{1, PS(kkk)}.upper_u];
                %         Delta_belief_add = Delta_belief_add*Partial_curve{kkk}.Belief_FE_belief(PS_add(kkk));
            end
        end

       
        bpa_2_remuve = 0;

       
        decomposition_end{num_sample}.step_two{num_combination}.bpa = bpa*bpa_scale - bpa_2_remuve;
%         
        final_FE_compact = [P_FE PS];
        final_FE = [position_u_uncoupled P_partial];
%         decomposition_end{num_sample}.step_two{num_combination}.FE = final_FE;
%                 decomposition_end{num_sample}.step_two{num_combination}.F = F - (num_functions_tot-1)*f_cost;
        decomposition_end{num_sample}.step_two{num_combination}.u_max = [upper_u_tot_decoupled coupled];%Partial_curve{1}.u_coupled(Sample{1, 1}.position(num_sample),:)];
%         
%         
%         
%         
        % table of FE reconstructed with this sample
        decomposition_end{num_sample}.final_FE = [decomposition_end{num_sample}.final_FE;...
            {final_FE_compact}...
            {final_FE}...
            F_belief(num_combination)...
            {decomposition_end{num_sample}.step_two{num_combination}.u_max}...
            decomposition_end{num_sample}.step_two{num_combination}.bpa...
            {lower_bound}...
            {upper_bound}];


    
    end
    
    
end

%--------------------------------------------------------------------------
% PLAUSIBILITY
%--------------------------------------------------------------------------

if in.output == 1 || in.output == 2
    
    d_plausibility = minmin.d;
    f_cost_Plausibility = problem_decomposition.objfun{n_obj}(d_plausibility, u_max_tot_Plausibility, problem_decomposition);
    
    
    % EVALUATE bpa-scale: for each sample a scaled curve is evaluated.
    
    bpa_scale_Plausibility = 1;
    for j = 1:length(PS)
        if in.dim_u_i(in.num_functions + j)>0
            if PS(j) == 1
                
                bpa_scale_Plausibility = bpa_scale_Plausibility*Sample{j}.Plausibility_f(PS(j));
            else
                
                bpa_scale_Plausibility = bpa_scale_Plausibility*(Sample{j}.Plausibility_f(PS(j)) - Sample{j}.Plausibility_f(PS(j)-1));
            end
        end
    end
    
    
    % RECONSTRUCTION of all the focal elements of the final scaled curves
    
    for num_combination = 1:num_combinations_max
        
        [P_FE] = position_sum_function(num_combinations_max, num_functions_tot, num_combination, decomposition_end, num_sample);
        
        bpa =1;
        
        F_Plausibility = 0;
        for i = 1:length(P_FE)
            bpa = bpa*decomposition_end{num_sample}.step_one{sub_not_empty(i)}{P_FE(i)}.bpa;
            
            F_Plausibility = F_Plausibility + decomposition_end{num_sample}.step_one{sub_not_empty(i)}{P_FE(i)}.downer_f;
        end
        
        
        
        decomposition_end{num_sample}.step_two{num_combination}.bpa_Plausibility = bpa*bpa_scale_Plausibility;
        
        
        
        decomposition_end{num_sample}.step_two{num_combination}.F_Plausibility = F_Plausibility - (num_functions_tot-1)*f_cost_Plausibility;
    end
    
    
end



end