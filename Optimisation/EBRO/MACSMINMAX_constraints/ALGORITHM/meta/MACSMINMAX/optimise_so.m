function [d,fval,exitflag,output] = optimise_so(problem_minmax, algo_outer, algo_inner, par_minmax)

global LR_inner


global GR_inner


global LR_outer


global GR_outer


global local_restart


global global_restart


global output_outer


global output_inner


global d_archive


global u_archive



% Rename inputs
n_d = problem_minmax.dim_d;
n_u = problem_minmax.dim_u;
n_obj = problem_minmax.n_obj;
sign_inner = problem_minmax.sign_inner; % 1 for minmax, -1 for minmin

nfevalmax = par_minmax.maxnfeval;



if sign_inner == -1

    problem_minmin = build_metaproblem_so_minmin(problem_minmax);
    par = algo_inner.par;
    par.nFeValMax = nfevalmax;     %% AD HOC FOR MPAIDEA!!!!!!!!!
    
    %% %%%%%%%%%%%%
    % CONSTRAINTS
    par.sign_inner=sign_inner;
    problem_minmin.par_objfun.sign = sign_inner;
    
    % objective and constraints are defined in different functions

    % Function to optimise
    fitnessfcn.obj       = problem_minmin.objfun;
    % Function of constraints
    if ~isempty(problem_minmin.par_objfun.constraint{1})
        fitnessfcn.constr    = problem_minmin.par_objfun.mask_constraints;
    else
        fitnessfcn.constr    = [];
    end
    
    
    % Flag to 0: objective and constraints are NOT in the same function
    fitnessfcn.obj_constr = 0;
    % How to handle constraints: set to 1 for weighted constraints with fixed
    % weights, or to 0 for penalty with no weights
    fitnessfcn.weighted = 0;
    % If the constraints are handled without weights, then define a tolerance
    % for the violation of the equality constraints
    fitnessfcn.ceq_eps = 1e-6;
    % Weights for penalty if fitnessfcn.weighted == 1
    fitnessfcn.w_ceq = 100;
    fitnessfcn.w_c = 100;

    
    problem_minmin.fitnessfcn = fitnessfcn;
    
    
    % Flag to 0: objective and constraints are NOT in the same function
%     problem_minmin.fitnessfcn.obj_constr = 0;
    % Weights for penalty
%     problem_minmin.fitnessfcn.w_ceq = 100;
%     problem_minmin.fitnessfcn.w_c   = 100;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%
    % FIXED PARAMETERS
    problem_minmin.par_objfun.problem_par_objfun{n_obj}.fix = problem_minmax.fix;
    %% %%%%%%%%%%%%
    
    [ dumin, fval_all_populations , exitflag , output_aux] = algo_inner.optimise(problem_minmin,par);
    
    [fval, min_pop] = min(fval_all_populations);
    dmin = dumin(min_pop, 1:n_d);
    umin = dumin(min_pop, n_d+1:end);
    d = dmin.*(problem_minmax.ub_d'-problem_minmax.lb_d') + problem_minmax.lb_d';
    u = cell(1,n_obj);
    obj=1;
    map_info = problem_minmin.par_objfun.map_u_info{obj};
    u{obj} = map_affine(umin,map_info);

    output.u = u;
    output.nfeval = output_aux.nfeval;
    % exitflag = 0;

else

    % local search flags, ADD SANITY CHECKS AND STUFF HERE
    lsflag_validation = par_minmax.local_search_flags.validation; % run either local search or func eval in validation
    lsflag_inner = par_minmax.local_search_flags.inner;           % run either local search or func eval in SO problem
    lsflag_outer = par_minmax.local_search_flags.outer;           % run either local search or func eval in MO problem

    %% %%%%%%%%%%%%
    % CONSTRAINTS i'm not sure it is useful
    algo_inner.par.sign_inner=sign_inner;
    algo_outer.par.sign_inner=sign_inner;
    %% %%%%%%%%%%%%

    
    % Build metaproblems
    problem_max_u = build_metaproblem_macsminmax_inner(problem_minmax);
    problem_min_d = build_metaproblem_macsminmax_outer(problem_max_u,lsflag_outer);


    %%%%%%%%%%%%%%
    % FIXED PARAMETERS
    problem_max_u.par_objfun.problem_par_objfun{n_obj}.fix=problem_minmax.fix;   
    %%%%%%%%%%%%%%
  
    
    % Random initial guess for d
    n_d0 = par_minmax.n_d0;
    d_0 = lhsu(zeros(1,n_d),ones(1,n_d),n_d0);

    % Initialise archive d
    d_record = [];
    f_record = [];
    
    f_constraint_record = [];
    % Now the fun begins
    nfeval = 0;

    
    % Initialise archive u    
    u_record = cell(1,n_obj);
    for i = 1:size(d_0,1)
        problem_max_u.par_objfun.d = d_0(i,:);                                                 % tell metaproblem what d to fix
        for obj = 1:n_obj
            problem_max_u.par_objfun.objective = obj;                                     % tell metaproblem what objective to optimise
            [ umax, ~ , ~ , output_aux ] = algo_inner.optimise(problem_max_u,algo_inner.par); % optimise
            u_record{obj} = [u_record{obj}; umax];                                        % append result
            nfeval = nfeval + output_aux.nfeval;
        end
    end

    
    
    
    
    
    
    
    %% MAIN LOOP
    
    size_u_record = n_obj;
    stop = false;
    while ~stop
        problem_min_d.par_objfun.u_record = u_record;
        
    %% %%%%%%%%%%%%
    % FIXED PARAMETERS
    problem_min_d.par_objfun.problem_fix_d.par_objfun.problem_par_objfun{n_obj}.fix = problem_minmax.fix;  
    problem_min_d.par_objfun.problem_par_objfun{n_obj}.fix=problem_minmax.fix;
    % ub_u for triangular problem 
    problem_min_d.par_objfun.problem_par_objfun{n_obj}.ub_u = problem_minmax.ub_u{1};
    %% %%%%%%%%%%%%        

    
        % "OUTER" LOOP: Compute dmin(i) = arg min {max f1(d,ue1), max f2(d,ue2)}      
        [dmin, f_outer, ~ , output_aux] = algo_outer.optimise(problem_min_d,algo_outer.par);
       
        GR_outer = [GR_outer output_aux.number_GR]; 
        LR_outer = [LR_outer output_aux.number_LR]; 
        
        output_outer = [output_outer; output_aux];

        nfeval = nfeval + (output_aux.nfeval+(n_obj>1))*(~lsflag_outer + lsflag_outer*20*problem_minmax.dim_u)*size_u_record; 
        % NOTE: the +1 is empirical
        % NOTE: 50*dim_u is hardcoded nfevalmax in local search but is too conservative and overestimates nfeval hence 20.
        
        dmin(dmin < 0) = 0;
        dmin(dmin > 1) = 1;

        
        % Remove solutions archived more than once (if any)
        [~,idx] = unique(round(1e8*dmin),'rows');
        dmin=dmin(idx,:);
        f_outer = f_outer(idx,:);
        f_record_aux = f_outer;

%         global fout
%         fout = [fout f_outer];
        
    
        
        % "INNER" LOOP: Compute ue(i+1) = arg max f(dmin(i),u)       
        for i = 1:size(dmin,1)
            problem_max_u.par_objfun.d = dmin(i,:);
            for obj = 1:n_obj
                
                % 1) evaluate the feasible  max(F(u))
                %
                
                problem_max_u.par_objfun.objective = obj;
                % Maximize subproblem to find umax = arg max f(dmin(i),u)
                [ umax, f_inner , ~ , output_aux] = algo_inner.optimise(problem_max_u,algo_inner.par);
                
        GR_inner = [GR_inner output_aux.number_GR]; 
        LR_inner = [LR_inner output_aux.number_LR]; 
 
        output_inner = [output_inner; output_aux];
                nfeval = nfeval + output_aux.nfeval;
                umax(umax < 0) = 0;
                umax(umax > 1) = 1;
                f_inner = -sign_inner*f_inner;
                
%         global fin
%         fin = [fin f_inner];


                %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % CONSTRAINTS
                if ~isempty(problem_max_u.fitnessfcn.constr)
                    
                    % 2) evaluate the max violation of the constraint
                    
%                     algo_inner.par;
%                     algo_inner.par.nFeValMax = 5000; % HC
                    [ umax_constraint, f_inner_constraint , ~ , output_aux] = optimise_constraint(problem_max_u,algo_inner.par);
                    % the sign  of the constraint function is changed inside
                    % "optimise_constraint": C is maximised minimizing - C
                    %
                    %               max(C) ---> min(-C)
                    
                    
                    nfeval = nfeval + output_aux.nfeval;
                    umax_constraint(umax_constraint < 0) = 0;
                    umax_constraint(umax_constraint > 1) = 1;
                    f_inner_constraint = sign_inner*f_inner_constraint;
                    
%                     problem_max_u.par_objfun.constraints_flag = 1;
                    %                 [F] = problem_max_u.objfun(u,problem_max_u.par_objfun);
                    
                    f_constraint_record = [f_constraint_record; f_inner_constraint];
                    % chose the unfeasible u if the constraint is not respected in all the domain.
                    epsilon = max(max(f_constraint_record)/3);
                    
                    % constraint relaxation
                    if any(f_inner_constraint > 0) && nfeval/nfevalmax < 0.7 ||...
                            any(f_inner_constraint > 0) + epsilon*(nfeval/nfevalmax > 0.7) + epsilon*(nfeval/nfevalmax > 0.8)
                        
                        f_inner = f_inner_constraint;
                        % if the constraint is violated i record that valoue
                        %
                        %
                        % the following is the function F valoue in the point u that
                        % maximises the constraint increased with a penalty
                        % function:
                        % - problem_max_u.objfun(umax_constraint, problem_max_u.par_objfun) + f_inner_constraint;
                        umax = [];
                        umax = umax_constraint;
                    end
                else
                    epsilon = 0;
                    f_inner_constraint = NaN;
                end
                %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                
                % Archive best between umax and umax2 [in case IDEA failed]
                if any(sign_inner*f_inner > sign_inner*f_outer(i,obj)) || any(f_inner_constraint > 0 + epsilon*(nfeval/nfevalmax > 0.6) + epsilon*(nfeval/nfevalmax > 0.8) + epsilon*(nfeval/nfevalmax > 0.9))
                    if lsflag_inner
                        u_cell_aux = cell(1,n_obj);
                        u_cell_aux{obj} = umax;
                        [fmax2,umax2,nfeval2] = u_validation(problem_max_u, dmin(i,:), u_cell_aux, lsflag_inner, obj);
                        nfeval = nfeval + nfeval2;

                        u_record{obj} = [u_record{obj}; umax2{obj}];
                        f_record_aux(i,obj) = fmax2(1,obj);
                    else
                        u_record{obj} = [u_record{obj}; umax];
                        clearvars f_record_aux
%                         f_record_aux(i,obj) = f_inner;
                        f_record_aux = f_inner;
                    end 

                    % If lsflag_inner = 0, the "else" step wastes evaluations and is
                    % useless: It finds and archives a u that is already in the
                    % archive, and is removed by the "unique" at line 96.
                else
                    if lsflag_inner
                        % Evaluate subproblem to find umax2 = arg max f(dmin(i),ua)
                        [fmax2,umax2,nfeval2] = u_validation(problem_max_u, dmin(i,:), u_record, lsflag_inner, obj);
                        nfeval = nfeval + nfeval2;
                        u_record{obj} = [u_record{obj}; umax2{obj}];
                        f_record_aux(i,obj) = fmax2(1,obj);
                    end
                end
            end
        end
        
        % Archive dmin
        d_record = [d_record; dmin];
        f_record = [f_record; f_record_aux];
        
        % Remove solutions archived more than once (if any)
        [~, idx] = unique(round(1e8*d_record),'rows');
        d_record = d_record(idx,:);
        f_record = f_record(idx,:); %not worrying too much by loss of validated solutions
        for obj = 1:n_obj
            [~, idx] = unique(round(1e8*u_record{obj}),'rows');
            u_record{obj} = u_record{obj}(idx,:);
        end
        
        size_u_record = sum(cellfun('size',u_record,1));
        nfeval_loop = nfevalmax - size_u_record*size(d_record,1)*(~lsflag_validation + lsflag_validation*20*problem_minmax.dim_u);
        % NOTE: 50*dim_u is hardcoded nfevalmax in local search but is too conservative and overestimates nfeval_val hence 20.

        
        stop = nfeval >= nfeval_loop;
    end

    %% archive cross check
    % % Refine solutions that maximize the objectives. This step also associates a u vector to each d vector
    % [ f_val, u_val_record, nfeval_val] = u_validation(problem_max_u, d_record, u_record, lsflag_validation, 1:n_obj);
    % nfeval = nfeval + nfeval_val;

    % %% Select non-dominated solutions and define outputs
    % sel = dominance(fval,0) == 0;

    % fval = f_val(sel,:);

    %% more clever archive cross check
    % even more clever would be to keep track of which u's have the d's been validated against already
    stop = false;
    checked = false(1,size(d_record,1));
    sel = false(1,size(d_record,1));
    u_val_record = cell(1,n_obj);
    for obj = 1:n_obj
        u_val_record{obj}=nan(size(d_record,1),problem_minmax.dim_u);
    end
    nfeval_val = 0;
    while ~stop
        sel = dominance(f_record,0) == 0;       % find non-dominated
        if(n_obj>1) 
            sel = sel';
        end
        
        tocheck = sel & ~checked;               % select only those that have not been checked
        stop = all(~tocheck);                   % stop if you have nothing to check
        % check those that have been selected
        d_tocheck = d_record(tocheck,:);
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % CONSTRAINTS
        [f_val_aux, u_val_record_aux, nfeval_aux] = u_validation_constraints(problem_max_u, d_tocheck, u_record, lsflag_validation, 1:n_obj);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        nfeval_val = nfeval_val + nfeval_aux;
        
        
        % update f_record and u_val_record
        f_record(tocheck,:) = f_val_aux;
        for obj = 1:n_obj
            u_val_record{obj}(tocheck,:) = u_val_record_aux{obj};
        end
        checked = checked | tocheck;           % update who has been checked
    end

    fval = f_record(sel,:);
    nfeval = nfeval + nfeval_val;

    %%
    frontsize = size(fval,1);

    d = d_record(sel,:).*repmat(problem_minmax.ub_d'-problem_minmax.lb_d',[frontsize,1]) + repmat(problem_minmax.lb_d',[frontsize,1]);
    u = cell(1,n_obj);
    for obj = 1:n_obj
        u_val_record{obj} = u_val_record{obj}(sel,:);
        map_info = problem_max_u.par_objfun.map_u_info{obj};
        for i = 1:frontsize
            u{obj}(i,:) = map_affine(u_val_record{obj}(i,:),map_info); %this can be easily vectorized
        end
    end

    output.u = u;
    output.nfeval = nfeval;
    exitflag = 0;
    
    
    d_archive = [d_archive; d_record];


    u_archive = [u_archive; u_record];
end




end