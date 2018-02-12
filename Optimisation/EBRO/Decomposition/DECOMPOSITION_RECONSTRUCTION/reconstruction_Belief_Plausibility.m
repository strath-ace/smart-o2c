function [decomposition_end, Plot_decomposition, num_sample_tot, LIST] = reconstruction_Belief_Plausibility(in, problem_decomposition, minmax, minmin, Sample, n_obj, n_point, algo_decomposition, Partial_curve)


%% here d and coupled u ARE FIXED
u_max_tot = 0;
u_max_tot_Plausibility = 0;


num_coupled_vectors = length(in.dim_u_i(in.num_functions +1 : end)); 


num_sample_tot = 1;

if in.output == 0 % evaluate only Belief
    
    d_belief = minmax.d;
    for index_num_coupled_vectors = 1:num_coupled_vectors
        if in.dim_u_i(in.num_functions + index_num_coupled_vectors)>0
            num_sample_tot = num_sample_tot*length(Sample{index_num_coupled_vectors}.FE);
        end
    end
    
    % input uncertain variable for Belief
    u_max_tot = minmax.u;
end

if in.output == 1  % evaluate only Plausibility
    
    d_plausibility = minmin.d;
    for index_num_coupled_vectors = 1:num_coupled_vectors
        if in.dim_u_i(in.num_functions + index_num_coupled_vectors)>0
            num_sample_tot = num_sample_tot*length(Sample{index_num_coupled_vectors}.FE_Plausibility);
        end
    end
    
    % input uncertain variable for Plausibility
    u_max_tot_Plausibility = minmin.u;
end

if in.output == 2  % evaluate both Belief and Plausibility
    
    d_belief = minmax.d;
    d_plausibility = minmin.d;
    
    for index_num_coupled_vectors = 1:num_coupled_vectors
        if in.dim_u_i(in.num_functions + index_num_coupled_vectors)>0
            num_sample_tot = num_sample_tot*length(Sample{index_num_coupled_vectors}.FE_Plausibility);
            if length(Sample{index_num_coupled_vectors}.FE) ~= length(Sample{index_num_coupled_vectors}.FE_Plausibility)
                disp('different number of sample between Belief and Plausibility')
            end
        end
    end
    % input uncertain variable for Belief and Plausibility
    u_max_tot = minmax.u;
    u_max_tot_Plausibility = minmin.u;
    
    
end

decomposition_end = cell(1,num_sample_tot);

%%


for num_sample = 1:num_sample_tot
    
    % vector of the positions of the sample in the h_ij-space (coupled)
    [PS] = position_Sample(num_sample_tot, num_coupled_vectors, num_sample, Sample, in); 
    
    
    for index_num_coupled_vectors = 1:num_coupled_vectors
        if in.dim_u_i(in.num_functions + index_num_coupled_vectors)>0
            
            if in.output == 0 || in.output == 2
                % BELIEF - components to optimize
                u_max_tot(var2opt(in.num_functions + index_num_coupled_vectors, problem_decomposition)) = Sample{index_num_coupled_vectors}.FE{PS(index_num_coupled_vectors)}.upper_u;
            end
            
            if in.output == 1 || in.output == 2
                % PLAUSIBILITY - components to optimize
                u_max_tot_Plausibility(var2opt(in.num_functions + index_num_coupled_vectors, problem_decomposition)) = Sample{index_num_coupled_vectors}.FE_Plausibility{PS(index_num_coupled_vectors)}.downer_u;
            end
        end
    end
    
    
    % problem_inner
    problem_inner = build_metaproblem_ideaminmax_s_inner(problem_decomposition);
    
    if in.output == 0 || in.output == 2
         problem_inner.par_objfun.d_belief = d_belief;
    end
    
    if in.output == 1 || in.output == 2
         problem_inner.par_objfun.d_plausibility = d_plausibility;
         
    end
    
    problem_inner.par_objfun.objective = n_obj;
    
    problem_inner.par_objfun.objfun = problem_decomposition.objfun;
    
    problem_inner.par_objfun.problem_par_objfun{n_obj} = problem_decomposition.par_objfun{n_obj};
    
    
    
    
    for i= 1:in.num_functions        % each u_i (uncoupled vector)
        
        if in.dim_u_i(i)>0
            
            
            
            problem_inner.par_objfun.vars_to_opt = var2opt(i,problem_decomposition);      
            problem_inner.dim = length(problem_inner.par_objfun.vars_to_opt);
            vars_to_opt = problem_inner.par_objfun.vars_to_opt;
            
            % init
            problem_max = problem_decomposition;
            
            problem_max.dim_u = problem_inner.dim;
            
            
            
            for index =1:problem_inner.dim
                lb_u(index,1) = {in.lb_u{n_obj}{vars_to_opt(index)}};
                ub_u(index,1) = {in.ub_u{n_obj}{vars_to_opt(index)}};
                bpa(index,1) = {in.bpa{n_obj}{vars_to_opt(index)}};
            end
            problem_max.lb_u = lb_u;
            problem_max.ub_u = ub_u;
            % BPA
            problem_max.bpa = bpa;
            
            
            
            
            
            
            problem_inner.bpa = problem_max.bpa;
            
            problem_decomposition.dim_u = problem_inner.dim;
            num_FE = number_Focal_Element(problem_decomposition, problem_inner);
            
            for ii=1:num_FE                                               % for each focal element
                problem_inner.lb = zeros(1,problem_inner.dim);
                problem_inner.ub = zeros(1,problem_inner.dim);
                
                position_FE = position(ii, problem_inner, problem_max);
                
                for iii = 1 : problem_inner.dim                            % find the domain of the fpcal element
                    
                    problem_inner.lb(iii) = problem_max.lb_u{iii,1}(position_FE(iii)); 
                    problem_inner.ub(iii) = problem_max.ub_u{iii,1}(position_FE(iii));
                    
                end
                
                
                [decomposition_end] = evaluate_max_min_FE_decoupled(i, ii, in, position_FE, decomposition_end, algo_decomposition, problem_inner, num_sample,u_max_tot, u_max_tot_Plausibility);
            end
            
        end
        
        
        
    end
    
    % max(F) = max(f_1)+max(f_2)+...
    [decomposition_end] = NEW_evaluate_combination_FE_decoupled(decomposition_end, num_sample, problem_decomposition, minmax, minmin, u_max_tot, Sample, PS, u_max_tot_Plausibility, n_obj, n_point, in, Partial_curve);
    
    
    
end


[LIST] = plot_Belief_Plausibility_decomposition(decomposition_end, in);



% if in.output == 0
%     % BELIEF
%     
%     
%     stairs(LIST.F_Bel, LIST.Bel,'linewidth',2)
%     xlabel('F')
%     ylabel('Belief ')
%     title('Final curve (decomposition approach)')
% end
% 
% 
% if in.output == 1
%     % PLAUSIBILITY
%      
%     figure
%     stairs(LIST.F_Pl, LIST.Pl,'linewidth',2)
%     xlabel('F')
%     ylabel('Plausibility')
%     title('Final curve (decomposition approach)')
% end
% 
% if in.output == 2
%     % BELIEF & PLAUSIBILITY
% 
%     figure
%     stairs(LIST.F_Bel, LIST.Bel,'linewidth',2)
%     hold on
%     stairs(LIST.F_Pl, LIST.Pl,'linewidth',2)
%     
%     xlabel('F')
%     ylabel('Belief & Plaiusibility')
%     title('Final curves (decomposition approach)')
% 
% end


for num_sample = 1:num_sample_tot
    if isempty(decomposition_end{num_sample})==0
        Plot_decomposition.max_f_min_f_FE{num_sample} = decomposition_end{num_sample}.step_two;
    end
end



end