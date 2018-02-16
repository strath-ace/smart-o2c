function problem = initialise_problem(dirname)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%
% This function creates the overall structure for multi-phase problems. It
% embeds all 'structures' from individual phases, creates an adjacency
% matrix and a cell array with all the linkage functions, so that it should
% be very easy to call the right linkage function even for arbitrarily
% complex multi-phase problems (not only sequential phases).

% Lorenzo A. Ricciardi, 2017

%% Create all structures

lb = [];
ub = [];
norm_lb = [];
norm_ub = [];

state_vars = logical([]);
control_vars = logical([]);
static_vars = logical([]);
other_vars = logical([]);
t0_vars = logical([]);
tf_vars = logical([]);
x0_vars = logical([]);
xf_vars = logical([]);
u0_vars = logical([]);
uf_vars = logical([]);
scales.xscale = [];
scales.uscale = [];
scales.tscale = [];
scales.other_vars_scale = [];
scales.scale_opt = [];

str = 'Initialising optimal control problem ...\n';
fprintf(str);

str = [dirname filesep 'phase_all.m'];
if ~(exist(str, 'file') == 2)
    
    error('No phases_all.m file is present');
    
end

run(str);

if ~exist('num_phases','var')
    
    error('num_phases not defined in phase_all.m');
    
else
    
    if num_phases<=0 || mod(num_phases,1)>0
        
        error('num_phases must be a positive integer');
        
    end
    
end

structure = cell(num_phases,1);
adj_matr = zeros(num_phases,num_phases);
link_fun_matr = cell(num_phases,num_phases);

if ~exist('constants','var')
    
    if num_phases>1
        
        warning('Multi-phase problem has no problem-level constants, but this is not explicitly stated. Include a constants=[] declaration in phases_all to suppress the warning');
        fprintf(str);
        
    end
    
    problem.constants = [];
    
else
    
    if num_phases == 1 && ~isempty(constants)
        
        error('Problem is single phase, but constants are defined at the problem level. This will be fine in the future, but for now include those constants in the phase definition.');
        
    else
        
        problem.constants = constants;
        
    end
    
end

if ~exist('other_vars_bounds','var')
    
    if num_phases>1
        
        warning('Multi-phase problem has no problem-level other_vars_bounds, but this is not explicitly stated. Include a other_vars_bounds=[] declaration in phases_all to suppress the warning');
        fprintf(str);
        
    end
    
    problem.other_vars_bounds = [];
    
else
    
    if num_phases == 1 && ~isempty(other_vars_bounds)
        
        error('Problem is single phase, but other_vars_bounds are defined at the problem level. This will be fine in the future, but for now include those constants in the phase definition.');
        
    else
        
        problem.other_vars_bounds = other_vars_bounds;
        
    end
        
end

if ~exist('other_vars_guess','var')
    
    if num_phases>1
        
        warning('Multi-phase problem has no problem-level other_vars_guess, but this is not explicitly stated. Include a other_vars_guess=[] declaration in phases_all to suppress the warning');
        fprintf(str);
        
    end
    
    problem.other_vars_guess = [];
    
else
    
    if num_phases == 1 && ~isempty(other_vars_guess)
        
        error('Problem is single phase, but other_vars_guess are defined at the problem level. This will be fine in the future, but for now include those constants in the phase definition.');
        
    else
        
        problem.other_vars_guess = other_vars_guess;
        
    end
        
end

if ~(size(problem.other_vars_guess,1)==size(problem.other_vars_bounds,1))
    
    error('problem.other_vars_guess must have the same number of rows as problem.other_vars_bounds');
    
end

if ~isempty(problem.other_vars_bounds)
    
    if any(problem.other_vars_guess<other_vars_bounds(:,1)) || any(problem.other_vars_guess>other_vars_bounds(:,2))
        
        error('Some values for other_vars_guess are outside other_vars_bounds');
        
    end
    
    problem.scales.other_vars_scale = max([problem.other_vars_bounds(:,2)-problem.other_vars_bounds(:,1),abs(problem.other_vars_bounds(:,1)),abs(problem.other_vars_bounds(:,2))],[],2);
    
end

if num_phases>1    
    
    if ~exist('g','var')
        
        error('No high level interphase objective functions defined');
        
    else
        
        problem.g = g;
        
    end
    
else
   
    if exist('g','var')
       
        error('Problem is defined as single phase. Objective function should be defined in phase_1 for better clarity. Remove it from phase_all');
        
    else
       
        problem.g = @single_phase_objectives;
                
    end
    
end

phase_mask = [];

for i = 1:num_phases
    
    str = ['Transcribing phase ', num2str(i) ' ...\n'];
    fprintf(str);
    
    str = [dirname filesep 'phase_',num2str(i)];
    
    if ~(exist(str, 'file') == 2)
        
        str2 = ['Could not locate ', str, '.m file'];
        error(str2);
        
    end
    
    run(str);
    
    %% Generation of transcription (problem independent!)
    
    structure{i} = prepare_transcription(num_eqs,num_controls,num_elems,state_order,control_order,DFET,state_distrib,control_distrib,test_distrib,integr_type);
    structure{i} = impose_boundary_conditions(structure{i},imposed_initial_states,imposed_final_states,imposed_t_0,imposed_t_f);
    
    %% Compute bounds, gauge normalisation factors
    
    structure{i} = transcribe_bounds(x_0,x_f,t_0,t_f,state_bounds,control_bounds,time_bounds,t0_bounds,tf_bounds,x0_bounds,xf_bounds,other_vars_bounds,other_vars_guess,control_guess,structure{i});
    
    %% Include function handles, derivatives, objective functions, constraints and weights
    
    structure{i} = include_functions(structure{i},f,dfx,dfu,g,weights,dgu0,dgxf,dguf,dgxi,dgui,c,dcx,dcu,e,dex,deu,h,wh,dhu0,dhxf,dhuf,dhxi,dhui,q,wq,dqu0,dqxf,dquf,dqxi,dqui,constants);
    
    %% Include structure constants, type of initialisation and info about next phases and the related linkage functions
    
    if ~(strcmp(init_type,'CloneIC') || strcmp(init_type,'DFET-all'))
            
            error('Unrecognized init_type: only ''CloneIC'' or ''DFET-all'' are implemented for now')
            
    else
    
        structure{i}.init_type = init_type;
        
    end
    
    
    
    structure{i}.next_phases = next_phases;
    
    structure{i}.constants = constants;
    
    if ~isempty(next_phases)
        
        structure{i}.link_funcs = link_funcs;
        
    end
    
    %% Stack all bounds in a single problem-wise vector
    
    lb = [lb; structure{i}.lbv];
    ub = [ub; structure{i}.ubv];
    norm_lb = [norm_lb; structure{i}.norm_lbv];
    norm_ub = [norm_ub; structure{i}.norm_ubv];
    
    %% Stack all flagging vectors, phase mask and scales in a single problem-wise vector
    
    state_vars = [state_vars; structure{i}.state_vars];
    control_vars = [control_vars; structure{i}.control_vars];
    static_vars = [static_vars; structure{i}.static_vars];
    other_vars = [other_vars; structure{i}.other_vars];
    t0_vars = [t0_vars; structure{i}.t0_vars];
    tf_vars = [tf_vars; structure{i}.tf_vars];
    x0_vars = [x0_vars; structure{i}.x0_vars];
    xf_vars = [xf_vars; structure{i}.xf_vars];
    u0_vars = [u0_vars; structure{i}.u0_vars];
    uf_vars = [uf_vars; structure{i}.uf_vars];
    
    phase_mask = [phase_mask; i*ones(size(structure{i}.lbv))];

    scales.xscale = [scales.xscale; structure{i}.scales.xscale];
    scales.uscale = [scales.uscale; structure{i}.scales.uscale];
    scales.tscale = [scales.tscale; structure{i}.scales.tscale];
    scales.other_vars_scale = [scales.other_vars_scale; structure{i}.scales.other_vars_scale];
    scales.scale_opt = [scales.scale_opt; structure{i}.scales.scale_opt];

end

%% Create adjacency matrix and matrix of linkage functions

for i = 1:num_phases
    
    nums = structure{i}.next_phases;
    
    % get if current phase is linked to another phase
    if ~isempty(nums)
        
        % check links make sense
        
        if any(nums<=0) || any(nums>num_phases) || any(mod(nums,1)>0)
            
            str = ['Values for next_phases in phase ',num2str(i),' are inconsistent. They must be positive integers, lower than the total nuber of phases (',num2str(num_phases),')'];
            
            error(str);
            
        else
            
            adj_matr(i,nums) = 1;
            
        end
        
        for j = 1:length(nums)
            
            if isempty(structure{i}.link_funcs{j})
                
                str = ['Missing linkage function between phase ',num2str(i), ' and phase ', num2str(j)];
                error(str);
                
            else
                
                link_fun_matr{i,nums(j)} = @(phase_a,phase_b) structure{i}.link_funcs{j}(phase_a,phase_b);
                
            end
            
        end
        
    end
    
end

%% Store structure(s), adjacency matrix and linkage functions in problem

problem.structure = structure;
problem.adj_matr = adj_matr;
problem.link_fun_matr = link_fun_matr;
problem.num_phases = num_phases;
problem.lb = lb;
problem.ub = ub;
problem.norm_lb = norm_lb;
problem.norm_ub = norm_ub;

problem.state_vars = state_vars;
problem.control_vars = control_vars;
problem.static_vars = static_vars;
problem.other_vars = other_vars;
problem.t0_vars = t0_vars;
problem.tf_vars = tf_vars;
problem.x0_vars = x0_vars;
problem.xf_vars = xf_vars;
problem.u0_vars = u0_vars;
problem.uf_vars = uf_vars;
problem.phase_mask = phase_mask;
problem.scales = scales;

end