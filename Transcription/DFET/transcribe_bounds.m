function structure = transcribe_bounds(x_0,x_f,t_0,t_f,state_bounds,control_bounds,time_bounds,t0_bounds,tf_bounds,x0_bounds,xf_bounds,other_vars_bounds,other_vars_guess,control_guess,structure)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%
% Automatically computes the scaling factors for all variables, given their
% bounds. Also structures the solution vector so that it always has the
% same ordering: generic static variables (if present), initial time 
% (if present), initial states (if present), initial controls (if present),
% all the normal variables, and final time (if present)

%% Input check

if size(state_bounds,1)~=structure.num_eqs
    
    error('state_bounds must have as many rows as equations')
    
else
    
    if size(state_bounds,2)~=2
        
        error('state_bounds must have 2 columns per row [lb ub]')
        
    else
        
        if any(state_bounds(:,1)>state_bounds(:,2))
            
            error('lower bounds for states must be on the first column of state_bounds')
            
        end
        
    end
    
end

if ~isempty(control_bounds)

if size(control_bounds,1)~=structure.num_controls
    
    error('control_bounds must have as many rows as controls')
    
else
    
    if size(control_bounds,2)~=2
        
        error('control_bounds must have 2 columns per row [lb ub]')
        
    else
        
        if any(control_bounds(:,1)>control_bounds(:,2))
            
            error('lower bounds for controls must be on the first column of control_bounds')
            
        end
        
    end 
    
end

if size(control_guess,2)>1
   
    error('control_guess must be a column vector');
    
else
    
    if size(control_guess,1)~=size(control_bounds,1)
    
        error('control_guess must have the same number of rows as control_bounds')
    
    else
        
        if any(control_guess<control_bounds(:,1)) || any(control_guess>control_bounds(:,2))
            
            error('some values of control_guess are out of bounds');
            
        end
        
    end        
    
end

end

if ~isequal(size(time_bounds), [1 2])
    
    error('time_bounds must be 1 row, 2 columns');
    
else
    
    if time_bounds(1)>=time_bounds(2)
        
        error('time_bounds(1) must be < than time_bounds(2)');
        
    end
    
end

if ~isempty(t0_bounds)
    
    if ~isequal(size(t0_bounds), [1 2])
        
        error('t0_bounds must be 1 row, 2 columns')
        
    else
        
        if t0_bounds(1)>t0_bounds(2)
            
            error('t0_bounds(1) must be greater than t0_bounds(2)')
            
        end
        
        if t0_bounds(1)<time_bounds(1)
            
            error('t0_bounds(1) must >= time_bounds(1)');
            
        end
        
        if t0_bounds(2)>time_bounds(2)
            
            error('t0_bounds(2) must >= time_bounds(1)');
            
        end
        
    end
    
else
    
    t0_bounds = time_bounds;
    
end

if ~isempty(tf_bounds)
    
    if ~isequal(size(tf_bounds), [1 2])
        
        error('tf_bounds must be 1 row, 2 columns')
        
    else
        
        if tf_bounds(1)>tf_bounds(2)
            
            error('tf_bounds(1) must be greater than tf_bounds(2)')
            
        end
        
        if tf_bounds(2)>time_bounds(2)
            
            error('tf_bounds(2) must be <= time_bounds(2)');
            
        end
        
    end
    
else
    
    tf_bounds = time_bounds;
    
end

if ~isempty(x0_bounds)
    
    if size(x0_bounds,1)>structure.num_eqs || size(x0_bounds,2)~=2
        
        error('x0_bounds must have at most num_eqs rows and exactly 2 columns')
        
    else
        
        if any(x0_bounds(:,1)>x0_bounds(:,2))
            
            error('lower bounds for x0 must be on the first column of x0_bounds')
            
        end
        
        if any(x0_bounds(:,1)<state_bounds(:,1))
            
            error('lower bounds for x0 must be >= lower state bounds');
            
        end
        
        if any(x0_bounds(:,2)>state_bounds(:,2))
            
            error('upper bounds for x0 must be <= lower state bounds');
            
        end
        
    end
    
else
    
    x0_bounds = state_bounds;
    
end

if ~isempty(xf_bounds)
    
    if size(xf_bounds,1)>structure.num_eqs || size(xf_bounds,2)~=2
        
        error('xf_bounds must have at most num_eqs rows and exactly 2 columns')
        
    else
        
        if any(xf_bounds(:,1)>xf_bounds(:,2))
            
            error('lower bounds for xf must be on the first column of xf_bounds')
            
        end
        
        if any(xf_bounds(:,1)<state_bounds(:,1))
            
            error('lower bounds for xf must be >= lower state bounds');
            
        end
        
        if any(xf_bounds(:,2)>state_bounds(:,2))
            
            error('upper bounds for xf must be <= lower state bounds');
            
        end
        
    end
    
else
    
    xf_bounds = state_bounds;
    
end

if ~isempty(other_vars_bounds)
    
    if size(other_vars_bounds,2)~=2
        
        error('other_vars_bounds must have 2 columns')
        
    else
        
        if any(other_vars_bounds(:,1)>other_vars_bounds(:,2))
            
            error('lower bounds for generic static variables must be on the first column of other_vars_bounds')
            
        end
        
    end    
    
end

if ~(size(other_vars_guess,1)==size(other_vars_bounds,1))
    
    error('other_vars_guess must have the same number of rows as other_vars_bounds');
    
end

%% Creation of ordered lower/upper bound vector for state and control vars

pad = zeros(structure.num_eqs*(structure.state_order+1)+structure.num_controls*(structure.control_order+1),2);
state_pad = pad(:,1);
control_pad = pad(:,1);

%% Creation of padding matrices, will simply be cloned for all elements

for i = 1:structure.num_eqs
    
    pad(1+(structure.state_order+1)*(i-1):(structure.state_order+1)*(i),:) = repmat(state_bounds(i,:),structure.state_order+1,1);
    state_pad(1+(structure.state_order+1)*(i-1):(structure.state_order+1)*(i),:) = 1;
    
end

for i = 1:structure.num_controls
    
    pad((structure.state_order+1)*structure.num_eqs+1+(structure.control_order+1)*(i-1):(structure.state_order+1)*structure.num_eqs+(structure.control_order+1)*(i),:) = repmat(control_bounds(i,:),structure.control_order+1,1);
    control_pad((structure.state_order+1)*structure.num_eqs+1+(structure.control_order+1)*(i-1):(structure.state_order+1)*structure.num_eqs+(structure.control_order+1)*(i),:) = 1;
    
end

%% Creation of lbv, ubv, and indexing vectors to identify which elements of the solution vector are state variables and (dynamic) control variables

lbv = repmat(pad(:,1),structure.num_elems,1);
ubv = repmat(pad(:,2),structure.num_elems,1);
state_vars = repmat(state_pad(:,1),structure.num_elems,1);
control_vars = repmat(control_pad(:,1),structure.num_elems,1);

% Bi-discontinuous DFET has additional variables (i.e. free final states)

if structure.DFET==1
    
    ids_free = 1:structure.num_eqs;
    ids_free = ids_free(structure.free_final_states);
    
    if ~isempty(ids_free)>0
        
        lbv = [lbv; xf_bounds(ids_free,1)];
        ubv = [ubv; xf_bounds(ids_free,2)];
        state_vars = [state_vars; ones(length(ids_free),1)];
        control_vars = [control_vars; zeros(length(ids_free),1)];
        
    end
        
end

%% Creation of flagging vectors for element ids, free xf, u0, uf

elem_id = repmat(1:structure.num_elems,size(pad,1),1);

elem_id = reshape(elem_id,size(elem_id,1)*size(elem_id,2),1);

% Bi-discontinuous DFET has additional variables (i.e. free final states)

if structure.DFET==1 && sum(structure.free_final_states)>0
       
    elem_id = [elem_id; structure.num_elems*ones(sum(structure.free_final_states),1)];
    
end

xf_vars = 0*lbv;
if sum(structure.free_final_states)>0
    
    xf_vars(end-sum(structure.free_final_states)+1:end) = 1;
    
end

u0_vars = (elem_id==1).*control_vars;

uf_vars = (elem_id==structure.num_elems).*control_vars;

%% Addition of t0, tf, x0 and generic static variables (if present)

addvars = [];
vars_type = [];

if ~isempty(other_vars_bounds)
    
    addvars = [addvars; other_vars_bounds];
    vars_type = [vars_type; ones(size(other_vars_bounds,1),1)];
    
    if any(other_vars_guess<other_vars_bounds(:,1)) || any(other_vars_guess>other_vars_bounds(:,2))
        
        error('Some values for other_vars_guess are outside other_vars_bounds');
        
    end
    
end

if ~structure.imposed_t0
    
    addvars = [addvars; t0_bounds];
    vars_type = [vars_type; 2];
    
    if isempty(t_0)
        
        error('No guess given for t_0');
        
    else
        
        if (t_0<t0_bounds(1)) || (t_0>t0_bounds(2))
                        
            error('Guess given for t_0 is out of bounds');
            
        end
        
    end
    
end

if any(~structure.imposed_initial_states)
    
    addvars = [addvars; x0_bounds(~structure.imposed_initial_states,:)];
    vars_type = [vars_type; 4*ones(sum(~structure.imposed_initial_states),1)];
    
    if any(x_0<x0_bounds(:,1)) || any(x_0>x0_bounds(:,2))
       
        error('Some initial guesses for x_0 are outside x0_bounds');
        
    end
end

if ~structure.imposed_tf
    
    addvars = [addvars; tf_bounds];
    vars_type = [vars_type; 3];
    
    if isempty(t_f)
        
        error('No guess given for t_f');
        
    else
        
        if (t_f<tf_bounds(1)) || (t_f>tf_bounds(2))
            
            error('Guess given for t_f is out of bounds');
            
        end
        
    end
    
end

if isempty(addvars)     % no static vars, thus no need to identify them
    
    t0_vars = 0*lbv;
    tf_vars = 0*lbv;
    x0_vars = 0*lbv;
    other_vars = 0*lbv;
    static_vars = 0*lbv;
    
else
    
    % trying to add all static variables at the beginning except for tf,
    % which would be at the end. This improves the sparsity pattern of the
    % Jacobian for multiphase problems
    
    lbv = [addvars(vars_type~=3,1); lbv; addvars(vars_type==3,1)];
    ubv = [addvars(vars_type~=3,2); ubv; addvars(vars_type==3,2)];
    
    old_vars = state_vars;
    state_vars = 0*lbv;
    state_vars(1+size(addvars(vars_type~=3,1)):end-size(addvars(vars_type==3,1))) = old_vars;
    
    old_vars = control_vars;
    control_vars = 0*lbv;
    control_vars(1+size(addvars(vars_type~=3,1)):end-size(addvars(vars_type==3,1))) = old_vars;
    
    static_vars = 0*lbv;
    static_vars(1:size(addvars(vars_type~=3,1))) = 1;
    
    if ~structure.imposed_tf
    
        static_vars(end) = 1;                                               % tf is a static variable
    
    end
    
    other_vars = 0*lbv;
    other_vars(vars_type==1) = 1;                                           %  miscellaneous static variables
    
    t0_vars = 0*lbv;
    t0_vars(vars_type==2) = 1;
    
    x0_vars = 0*lbv;
    x0_vars(vars_type==4) = 1;
    
    tf_vars = 0*lbv;
    if ~structure.imposed_tf
    
        tf_vars(end) = 1;
    
    end 
    
    old_vars = elem_id;
    elem_id = 0*lbv;
    elem_id(1+size(addvars(vars_type~=3,1)):end-size(addvars(vars_type==3,1))) = old_vars;

    old_vars = xf_vars;
    xf_vars = 0*lbv;
    xf_vars(1+size(addvars(vars_type~=3,1)):end-size(addvars(vars_type==3,1))) = old_vars;
        
    old_vars = u0_vars;
    u0_vars = 0*lbv;
    u0_vars(1+size(addvars(vars_type~=3,1)):end-size(addvars(vars_type==3,1))) = old_vars;
    
    old_vars = uf_vars;
    uf_vars = 0*lbv;
    uf_vars(1+size(addvars(vars_type~=3,1)):end-size(addvars(vars_type==3,1))) = old_vars;
     
end

% convert to logical arrays
state_vars = logical(state_vars);
control_vars = logical(control_vars);
static_vars = logical(static_vars);
t0_vars = logical(t0_vars);
tf_vars = logical(tf_vars);
x0_vars = logical(x0_vars);
other_vars = logical(other_vars);
xf_vars = logical(xf_vars);
u0_vars = logical(u0_vars);
uf_vars = logical(uf_vars);

%% Computation of normalisation factors

% states, controls and time scales

xscale = max([state_bounds(:,2)-state_bounds(:,1),abs(state_bounds(:,1)),abs(state_bounds(:,2))],[],2);

if ~isempty(control_bounds)

    uscale = max([control_bounds(:,2)-control_bounds(:,1),abs(control_bounds(:,1)),abs(control_bounds(:,2))],[],2);
else
    
    uscale = [];
    
end

tscale = max([time_bounds(:,2)-time_bounds(:,1),abs(time_bounds(:,1)),abs(time_bounds(:,2))],[],2);

if ~isempty(other_vars_bounds)
    
    other_vars_scale = max([other_vars_bounds(:,2)-other_vars_bounds(:,1),abs(other_vars_bounds(:,1)),abs(other_vars_bounds(:,2))],[],2);
    
else
    
    other_vars_scale = zeros(0,1);
    
end

% optimisation vector scale

scale_opt = max([ubv-lbv,abs(ubv),abs(lbv)],[],2); % this should ensure that scale_opt is always >= lbv and <= ubv and avoid normalisation errors

% correction of optimisation vector scale for final time variables (needed
% if final time bounds are narrower than state bounds)

if any(structure.free_final_states) && (structure.DFET==1)
    
    nfree = sum(structure.free_final_states);
    
    scale_opt(end-nfree+1-any(vars_type==3):end-any(vars_type==3)) = max([scale_opt(end-nfree+1-any(vars_type==3):end-any(vars_type==3)),xscale(structure.free_final_states)],[],2);
    %scale_opt(end-nfree+1:end) = max([scale_opt(end-nfree+1:end),xscale(structure.free_final_states)],[],2);
        
end


% correction of optimisation vector scale for static variables (boundary
% times and states)

if ~isempty(scale_opt(t0_vars))
    
    scale_opt(t0_vars) = max(scale_opt(t0_vars),tscale);
    
end

if ~isempty(scale_opt(tf_vars))
    
    scale_opt(tf_vars) = max(scale_opt(tf_vars),tscale);
    
end

if ~isempty(scale_opt(x0_vars))
    
    scale_opt(x0_vars) = max([scale_opt(x0_vars),xscale(~structure.imposed_initial_states)] , [] , 2);
    
end

norm_lbv = lbv./scale_opt;
norm_ubv = ubv./scale_opt;

%% Addition of all bounds, scaling factors and variable identity flags to transcription

structure.state_bounds = state_bounds;
structure.control_bounds = control_bounds;
structure.time_bounds = time_bounds;

if ~isempty(x0_bounds)
    
    structure.x0_bounds = x0_bounds;
    
else
    
    structure.x0_bounds = state_bounds;
    
end

if ~isempty(xf_bounds)
    
    structure.xf_bounds = xf_bounds;
    
else
    
    structure.xf_bounds = state_bounds;
    
end

if ~isempty(t0_bounds)
    
    structure.t0_bounds = t0_bounds;
    
else
    
    structure.t0_bounds = time_bounds;
    
end

if ~isempty(tf_bounds)
    
    structure.tf_bounds = tf_bounds;
    
else
    
    structure.tf_bounds = time_bounds;
    
end

if ~isempty(other_vars_bounds)
    
    structure.other_vars_bounds = other_vars_bounds;
    
else
    
    structure.other_vars_bounds = [];
    
end

structure.other_vars_bounds = other_vars_bounds;
structure.lbv = lbv;
structure.ubv = ubv;
structure.scales.xscale = xscale;
structure.scales.uscale = uscale;
structure.scales.tscale = tscale;
structure.scales.other_vars_scale = other_vars_scale;
structure.scales.scale_opt = scale_opt;
structure.norm_lbv = norm_lbv;
structure.norm_ubv = norm_ubv;
structure.state_vars = state_vars;
structure.control_vars = control_vars;
structure.static_vars = static_vars;
structure.t0_vars = t0_vars;
structure.tf_vars = tf_vars;
structure.x0_vars = x0_vars;
structure.xf_vars = xf_vars;
structure.other_vars = other_vars;
structure.t_0 = t_0;
structure.t_0_norm = t_0/tscale;
structure.t_f = t_f;
structure.t_f_norm = t_f/tscale;
structure.x_0 = x_0;
structure.x_0_norm = x_0./xscale;
structure.x_f = x_f;
structure.x_f_norm = x_f./xscale;
structure.other_vars_guess = other_vars_guess;
structure.control_guess = control_guess;
structure.elem_id = elem_id;
structure.u0_vars = u0_vars;
structure.uf_vars = uf_vars;

end