function [algo_minmin] = map_u_info_minmin(minmax_problem, obj, problem_minmin_cost)

% parameters

%     if isfield(minmax_problem,'maxnfeval')
%         par_minmax.maxnfeval = minmax_problem.maxnfeval; % function evaluation limit (does not account for ls in validation)
%     else
%         error('max nfeval not supplied')
%     end

    par_minmax.local_search_flags.validation = false;
    par_minmax.keep_d_record = true;
    par_minmax.tol_conv = 1e-3; %Tolerance to convergence, only used in S.O. case
    par_minmax.indicator_d_max = 1e-3; % Threshold on the indicator (EI, PI, etc.) for outer loop
    par_minmax.indicator_u_max = 1e-3; % Threshold on the indicator (EI, PI, etc.) for inner loop
    par_minmax.iter_d_max = 20; %20*minmax_problem.dim_d; % Max. iterations for outer loop in an iteration of MinMaReK
    par_minmax.iter_u_max = 20; %20*minmax_problem.dim_u; % Max. iterations for inner loop in an iteration of MinMaReK

    
    % initialisation of the surrogate parameters
        if minmax_problem.n_obj == 1
            par_minmax.indicator_d = str2func('kriging_EI2');
        elseif minmax_problem.n_obj == 2
            par_minmax.indicator_d = str2func('kriging_EIaug');
        else
            error('We do not have kriging indicators for many-objective (yet)' );
        end
        
        info_minmin = cell(1,minmax_problem.n_obj);
%         for obj = 1:minmax_problem.n_obj
            map_u_info = get_map_info(minmax_problem.lb_u{obj}, minmax_problem.ub_u{obj});

            num_fe = 1;
            for dim = 1:problem_minmin_cost.dim_u
                num_fe = num_fe * map_u_info.n_int{dim};
            end

            info_minmin{obj} = struct;
            info_minmin{obj}.dim_d = minmax_problem.dim_d;
            info_minmin{obj}.dim_u = minmax_problem.dim_u;
            info_minmin{obj}.method = 'kriging';
            info_minmin{obj}.corrfun = @corrgauss;
            info_minmin{obj}.regrfun = @regpoly0;
            info_minmin{obj}.training = str2func([lower(info_minmin{obj}.method) '_training']);
            info_minmin{obj}.predictor = str2func([lower(info_minmin{obj}.method) '_predictor']);

            info_minmin{obj}.indicator_u = str2func([lower(info_minmin{obj}.method) '_EI']);
            info_minmin{obj}.model = cell(1,num_fe);
            info_minmin{obj}.changed_FE = zeros(1,num_fe);
            info_minmin{obj}.x_doe = cell(1,num_fe);
            info_minmin{obj}.f_doe = cell(1,num_fe);
            info_minmin{obj}.map_info = map_u_info;
            
            info_minmin{obj}.find = @surrogate_find_fe;
            info_minmin{obj}.add = @surrogate_add;
            info_minmin{obj}.update = @surrogate_update;
            info_minmin{obj}.num_FE = num_fe;

%         end
    par_minmax.surrogate = info_minmin;

algo_minmin.par_minmax = par_minmax;

end