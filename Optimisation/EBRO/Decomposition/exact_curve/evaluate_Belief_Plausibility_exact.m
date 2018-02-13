function [EXACT_FE, LIST_EXACT] = evaluate_Belief_Plausibility_exact(in, problem_decomposition_TEST, minmax, minmin, n_obj, algo_decomposition)


%% freeze variables u and run maximization
problem_inner = build_metaproblem_ideaminmax_s_inner(problem_decomposition_TEST);

if in.output == 0 || in.output == 2
    problem_inner.par_objfun.d_belief = minmax.d;
    problem_inner.par_objfun.u_belief = minmax.u;
end

if in.output == 1 || in.output == 2
    problem_inner.par_objfun.d_plausibility = minmin.d;
    problem_inner.par_objfun.u_plausibility = minmin.u;
end






% problem_inner.par_objfun.d = dmin;
problem_inner.par_objfun.objective = n_obj;
% problem_inner.par_objfun.umax = u_max_tot;
% problem_inner.objfun = @mask_objfun_max_decomposition;
problem_inner.par_objfun.objfun = problem_decomposition_TEST.objfun;

problem_inner.par_objfun.problem_par_objfun{n_obj} = problem_decomposition_TEST.par_objfun{n_obj};



% structure decomposition
num_FE = number_Focal_Element_TEST(problem_decomposition_TEST);
EXACT_FE = cell(1,num_FE);




problem_inner.par_objfun.vars_to_opt = [1:problem_inner.dim];              % all the components of vector u
problem_inner.dim = length(problem_inner.par_objfun.vars_to_opt);

% init
        problem_max = problem_decomposition_TEST;
        
        problem_max.dim_u = problem_inner.dim;

        
        for index =1:problem_inner.dim
                lb_u(index,1) = {in.lb_u{n_obj}{index}};
                ub_u(index,1) = {in.ub_u{n_obj}{index}};
                bpa(index,1) = {in.bpa{n_obj}{index}};
        end
        problem_max.lb_u = lb_u;
        problem_max.ub_u = ub_u;
        % BPA
        problem_max.bpa = bpa;
        



% initialize algorithm
problem_inner.n_obj = n_obj;
problem_inner.dim_d = in.dim_d;
problem_inner.dim_u = sum(in.dim_u);


problem_inner.bpa = problem_max.bpa;

problem_decomposition_TEST.dim_u = problem_inner.dim;


LIST_EXACT.FE = [];
for ii=1:num_FE
    problem_inner.lb = zeros(1,problem_inner.dim);
    problem_inner.ub = zeros(1,problem_inner.dim);
    
    position_FE = position(ii, problem_inner, problem_max);
    
    for iii = 1 : problem_inner.dim
        
        problem_inner.lb(iii) = problem_max.lb_u{iii,1}(position_FE(iii));
        problem_inner.ub(iii) = problem_max.ub_u{iii,1}(position_FE(iii));
        
    end
    
    [EXACT_FE] = max_func_decomposizion_TEST(in, ii, position_FE, EXACT_FE, algo_decomposition, problem_inner);
    
    LIST_EXACT.FE = [LIST_EXACT.FE;...
                     {position_FE}...
                     EXACT_FE{ii}.upper_f...
                     EXACT_FE{ii}.upper_u...
                     EXACT_FE{ii}.bpa...
                     {problem_inner.lb}...
                     {problem_inner.ub}];
                 
                 
%             savefolder = strcat('RESULTS_spacecraft_1_16000/exact/');
%             mkdir(savefolder);
%             save(strcat(savefolder,'spacecraft_exact_1_16000'));


            
end


global num_maximization_exact;
num_maximization_exact = num_FE;


%% plot Bl

for k = 1 : length(EXACT_FE)     % vector of all the f-values of the FEs
    if in.output == 0 || in.output == 2
       f_belief(k) = EXACT_FE{k}.upper_f;
    end
    
    if in.output == 1 || in.output == 2
       f_plausibility(k) = EXACT_FE{k}.downer_f;
    end
    
    bpa_end(k) = EXACT_FE{k}.bpa;
end


if in.output == 0 || in.output == 2
    
    [f_sorted_Bel, position_sorted_Bel] = sort(f_belief);
    
    LIST_EXACT.F_Bel = f_sorted_Bel;
    LIST_EXACT.Bel   = cumsum(bpa_end(position_sorted_Bel));
end


% Plausibility
if in.output == 1 || in.output == 2
    
    [f_sorted_Pl, position_sorted_Pl] = sort(f_plausibility);
    
    LIST_EXACT.F_Pl = f_sorted_Pl;
    LIST_EXACT.Pl   = cumsum(bpa_end(position_sorted_Pl));
end




end