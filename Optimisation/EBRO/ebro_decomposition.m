% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [minmin, minmax, LIST, LIST_EXACT] = ebro_decomposition(problem, algo_minmax, algo_outer, algo_inner, algo_decomposition)
%
% 


minmax = [];
minmin = [];
LIST = [];
LIST_EXACT = [];



if problem.input == 0
    
    if problem.output == 0 || problem.output == 2
        [minmax] = evaluate_minmax(problem, algo_minmax, algo_outer, algo_inner);
        n_points_tot = length(minmax.d(:,1));
    end
    
    if problem.output == 1 || problem.output == 2
        [minmin] = evaluate_minmin(problem, algo_minmax, algo_outer, algo_inner);
        n_points_tot = length(minmin.d(:,1));
    end
    
    
    
    
    n_obj_tot = problem.n_obj;
    
    minmax_input = minmax;
    minmin_input = minmin;
    problem_input = problem;
    
    
elseif problem.input == 1
    
    minmax = [];
    minmin = [];
    
    if isempty(problem.input_minmin_minmax)
        
        minmax_input.d = problem.input_d_max;
        minmax_input.u = problem.input_u_max;
        minmax_input.f = problem.input_F_max;
        
        minmin_input.d = problem.input_d_min;
        minmin_input.u = problem.input_u_min;
        minmin_input.f = problem.input_F_min;
        
    else
        
        load(problem.input_minmin_minmax, 'minmax', 'minmin')
        
        minmax_input = minmax;
        minmin_input = minmin;
    end
    
    
    problem_input = problem;
    n_points_tot = length(minmax.d(:,1));
    n_obj_tot = problem.n_obj;
    
elseif problem.input == 2                                     % input d, run u
    
    minmax = [];  
    minmin = []; 
    max    = []; 
    min    = [];
    
    
    if isempty(problem.input_minmin_minmax)
        if problem.output == 0 || problem.output == 2
            minmax.d = problem.input_d_max;
        end
        if problem.output == 1 || problem.output == 2
            minmin.d = problem.input_d_min;
        end
    else
        load(problem.input_minmin_minmax, 'minmax', 'minmin')
    end
    
    
    if problem.output == 0 || problem.output == 2
        problem.sign_inner = 1;
        [max] = evaluate_max(problem, minmax.d, algo_inner);
    end
    
    if problem.output == 1 || problem.output == 2
        problem.sign_inner = -1;
        [min] = evaluate_min(problem, minmin.d, algo_inner);
    end
    
    n_points_tot = 1;
    n_obj_tot = 1;
    
    minmax_input = max;
    minmin_input = min;
    problem_input = problem;
    
end



for n_point = 1:n_points_tot                  % all the design points(DP); if single objective there is only one DP.
    
    for n_obj = 1:n_obj_tot
        
        if isfield(minmax_input,'d')
            minmax.d = minmax_input.d(n_point,:);
            minmax.u = minmax_input.u{n_obj}(n_point,:);
            minmax.f = minmax_input.f(n_point, n_obj);
        end
        if isfield(minmin_input,'d')
            minmin.d = minmin_input.d(n_point,:);
            minmin.u = minmin_input.u{n_obj}(n_point,:);
            minmin.f = minmin_input.f(n_point, n_obj);
        end
        
        %% DECOMPOSITION
        problem.sign_inner = 1;
        
        if problem.exact_curves ~=2
            % PARTIAL CURVES (coupled variables)
            [decomposition, Partial_curve] = evaluate_Belief_Plausibility_coupled_vectors(problem, problem, minmax, minmin, n_obj, n_point, algo_decomposition);
            
            % SAMPLE
            [Sample] = Sampling_Belief_Plausibility(Partial_curve, decomposition, problem);
            
            % RECONSTRUCTION
            [decomposition_end, Plot_decomposition, num_sample_tot, LIST] = reconstruction_Belief_Plausibility(problem, problem, minmax, minmin, Sample, n_obj, n_point, algo_decomposition, Partial_curve);
            
            
            if problem.exact_curves == 0
                hold on
                stairs(LIST.F_Bel, LIST.Bel,'b', 'linewidth',2)
            end
        end
        
        
        
        %% EXACT CURVES (Belief and/or Plausibility)
        if problem.exact_curves == 1
            
            [EXACT_FE, LIST_EXACT] = evaluate_Belief_Plausibility_exact(problem, problem, minmax, minmin, n_obj, algo_decomposition);
            
            
            hold on
            if problem.output == 0
                stairs(LIST_EXACT.F_Bel, LIST_EXACT.Bel, 'k', 'linewidth', 1)
                stairs(LIST.F_Bel, LIST.Bel,'b', 'linewidth',2)
                legend('Belief Exact','Belief Decomposition')
            elseif problem.output == 1
                stairs(LIST_EXACT.F_Pl, LIST_EXACT.Pl, 'k', 'linewidth',1)
                stairs(LIST.F_Pl, LIST.Pl,'r','linewidth',2)
                legend('Plausibility Exact','Plausibility Decomposition')
            elseif problem.output == 2
                stairs(LIST_EXACT.F_Bel, LIST_EXACT.Bel, 'k', 'linewidth', 1)
                stairs(LIST_EXACT.F_Pl, LIST_EXACT.Pl, 'k', 'linewidth', 1)
                stairs(LIST.F_Bel, LIST.Bel,'b','linewidth',2)
                stairs(LIST.F_Pl, LIST.Pl,'r','linewidth',2)
                legend('Belief e Plausibility Exact','Belief e Plausibility Decomposition')
            end
            hold off
            
            
        elseif problem.exact_curves == 2
            [EXACT_FE, LIST_EXACT] = evaluate_Belief_Plausibility_exact(problem, problem, minmax, minmin, n_obj, algo_decomposition);
            
            
            hold on
            if problem.output == 0
                stairs(LIST_EXACT.F_Bel, LIST_EXACT.Bel, 'b', 'linewidth', 3)
                legend('Belief Exact','Belief Decomposition')
            elseif problem.output == 1
                stairs(LIST_EXACT.F_Pl, LIST_EXACT.Pl, 'r', 'linewidth', 3)
                legend('Plausibility Exact','Plausibility Decomposition')
            elseif problem.output == 2
                stairs(LIST_EXACT.F_Bel, LIST_EXACT.Bel, 'b', 'linewidth', 3)
                stairs(LIST_EXACT.F_Pl, LIST_EXACT.Pl, 'r', 'linewidth', 3)
                legend('Belief e Plausibility Exact','Belief e Plausibility Decomposition')
            end
            hold off
            
        end
        
        
        
        
    end
    
end


return