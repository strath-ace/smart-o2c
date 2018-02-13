function [minmin, minmax, LIST, LIST_EXACT] = ebro_decomposition(problem, algo_minmax, algo_outer, algo_inner, algo_decomposition)



if problem.input == 0
    
minmax = [];
minmin = [];
    
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
    
    load('minmax_minmin_TC1', 'minmax', 'minmin')
    
    n_points_tot = length(minmax.d(:,1));
    n_obj_tot = problem.n_obj;
    
    minmax_input = minmax;
    minmin_input = minmin;
    problem_input = problem;
    
elseif problem.input == 2                                     % input d, run u
    %     problem.maxnfeval = 1e5;
    
    minmax = [];  minmin = []; max = []; min = [];
    
    
    %     load('minmax_TC1_ok')
    load('minmax_minmin_TC1', 'minmax', 'minmin')
    d=minmax.d;
    
    if problem.output == 0 || problem.output == 2
        [max] = evaluate_max(problem, d, algo_inner);
    end
    
    if problem.output == 1 || problem.output == 2
        [min] = evaluate_min(problem, d, algo_inner);
    end
    
    n_points_tot = 1;
    n_obj_tot = 1;
    
    minmax_input = max;
    minmin_input = min;
    problem_input = problem;
    
end



for n_point = 1:n_points_tot                  % all the design points
    
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
        
 problem.sign_inner = 1; 
        %% DECOMPOSITION
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
 
 
        
        %% EXACT BELIEF CURVE
        if problem.exact_curves == 1      
            [EXACT_FE, LIST_EXACT] = evaluate_Belief_Plausibility_exact(problem, problem, minmax, minmin, n_obj, algo_decomposition);
            
            
            hold on
            if problem.output == 0
                stairs(LIST_EXACT.F_Bel, LIST_EXACT.Bel, 'b', 'linewidth', 3)
                stairs(LIST.F_Bel, LIST.Bel,'b', 'linewidth',2)
                legend('Belief Exact','Belief Decomposition')
            elseif problem.output == 1
                stairs(LIST_EXACT.F_Pl, LIST_EXACT.Pl, 'r', 'linewidth', 3)
                stairs(LIST.F_Pl, LIST.Pl,'r','linewidth',2)
                legend('Plausibility Exact','Plausibility Decomposition')
            elseif problem.output == 2
                stairs(LIST_EXACT.F_Bel, LIST_EXACT.Bel, 'b', 'linewidth', 3)
                stairs(LIST_EXACT.F_Pl, LIST_EXACT.Pl, 'r', 'linewidth', 3)                
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