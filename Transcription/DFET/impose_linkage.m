function [c,Jc] =  impose_linkage(phase_boundaries,this_phase_id,next_phase_id,problem,jacflag)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%
% This function computes all linkage conditions of this phase, and their
% Jacobian. 

%% Compute constraints

% create phases structures


% compute generic linkage function

this_c = problem.link_fun_matr{this_phase_id,next_phase_id}(phase_boundaries{this_phase_id},phase_boundaries{next_phase_id});

% apopend all linkage constraints (time matching at the end, since time is
% the last variable of each phase and it's tidier to do it this way)

c = [this_c];

%% Compute all Jacobians

% initialise Jacobians as empty matrices of the correct "size" (avoids
% problems at the outer loop) 

Jc.static_a = zeros(0,size(c,1));
Jc.x0_a = zeros(0,size(c,1));
Jc.u0_a = zeros(0,size(c,1));
Jc.t0_a = zeros(0,size(c,1));
Jc.xf_a = zeros(0,size(c,1));
Jc.uf_a = zeros(0,size(c,1));
Jc.tf_a = zeros(0,size(c,1));

Jc.static_b = zeros(0,size(c,1));
Jc.x0_b = zeros(0,size(c,1));
Jc.u0_b = zeros(0,size(c,1));
Jc.t0_b = zeros(0,size(c,1));
Jc.xf_b = zeros(0,size(c,1));
Jc.uf_b = zeros(0,size(c,1));
Jc.tf_b = zeros(0,size(c,1));


if jacflag == 1
    
    h = 1e-7;
    
    %% Derivatives wrt phase A variables
    
    phase_sel = this_phase_id;
    
    % derivatives wrt static variables (phase A)
    
    nstatic_a = length(phase_boundaries{phase_sel}.static);                  % get number of static variables of this phase
    
    if nstatic_a>0
        
        Jstatic = zeros(size(c,1),nstatic_a);
        
        for i=1:nstatic_a
            
            temp_static = phase_boundaries{phase_sel}.static;               % get static variables of this phase
            temp_static(i) = temp_static(i)+h;                              % add small increment
            
            phase_boundaries_temp = phase_boundaries;                       % clone phase_boundaries structure
            phase_boundaries_temp{phase_sel}.static = temp_static;          % modify temporary clone, adding the incremented static variables
            
            c2 = impose_linkage(phase_boundaries_temp,this_phase_id,next_phase_id,problem,0);   % call this function, recursively, to get the values
            Jstatic(:,i) = (c2-c)'/h;                                       % compute finite difference approx of Jacobian
            
        end
        
        Jc.static_a = Jstatic';                                             % fmincon uses the Transposed of the Jacobian...
        
    end
    
    % derivatives wrt t0 (phase A)
    
    if problem.structure{phase_sel}.imposed_t0 == 0
        
        temp_t0 = phase_boundaries{phase_sel}.t0;                           % get t0 variables of this phase
        temp_t0 = temp_t0+h;                                                % add small increment
        
        phase_boundaries_temp = phase_boundaries;                           % clone phase_boundaries structure
        phase_boundaries_temp{phase_sel}.t0 = temp_t0;                      % modify temporary clone, adding the incremented static variables
        
        c2 = impose_linkage(phase_boundaries_temp,this_phase_id,next_phase_id,problem,0);   % call this function, recursively, to get the values
        Jt0 = (c2-c)/h;                                                     % compute finite difference approx of Jacobian
        
        Jc.t0_a = Jt0';
        
    end
    
    % derivatives wrt x0 (phase A)
    
    ids = 1:length(phase_boundaries{phase_sel}.x0);                         % indexing vector for id of free x0
    ids = ids(~problem.structure{phase_sel}.imposed_initial_states);        % keep only indexes of free x0
    nx0 = length(ids);                                                      % number of free initial states of this phase
    
    if nx0>0
        
        Jx0 = zeros(size(c,1),nx0);
        
        for i=1:nx0
            
            temp_x0 = phase_boundaries{phase_sel}.x0;                       % get initial state variables of this phase
            temp_x0(ids(i)) = temp_x0(ids(i))+h;                            % add small increment
            
            phase_boundaries_temp = phase_boundaries;                       % clone phase_boundaries structure
            phase_boundaries_temp{phase_sel}.x0 = temp_x0;                  % modify temporary clone, adding the incremented static variables
            
            c2 = impose_linkage(phase_boundaries_temp,this_phase_id,next_phase_id,problem,0);   % call this function, recursively, to get the values
            Jx0(:,i) = (c2-c)'/h;                                            % compute finite difference approx of Jacobian
            
        end
        
        Jc.x0_a = Jx0';
        
    end
    
    % derivatives wrt u0 (phase A)
    
    nu0 = length(phase_boundaries{phase_sel}.u0);                           % get number of control variables in this phase
    
    if nu0>0
        
        Ju0 = zeros(size(c,1),nu0);
        
        for i=1:nu0
            
            temp_u0 = phase_boundaries{phase_sel}.u0;                       % get final control VALUES of this phase
            temp_u0(i) = temp_u0(i)+h;                                      % add small increment
            
            phase_boundaries_temp = phase_boundaries;                       % clone phase_boundaries structure
            phase_boundaries_temp{phase_sel}.u0 = temp_u0;                  % modify temporary clone, adding the incremented static variables
            
            c2 = impose_linkage(phase_boundaries_temp,this_phase_id,next_phase_id,problem,0);   % call this function, recursively, to get the values
            Ju0(:,i) = (c2-c)'/h;                                        % compute finite difference approx of Jacobian
            
        end
        
        % SUPER IMPORTANT: NEED TO CAST THIS DERIVATIVE IN THE SPACE OF
        % THE ACTUAL VARIABLES, SO WE NEED TO MULTIPLY THOSE
        % DERIVATIVES BY THE BASIS FUNCTIONS OF THE CONTROLS (LIKE FOR
        % THE DERIVATIVES OF THE OTHER CONSTRAINTS)
        
        Jc.u0_a = (Ju0*problem.structure{phase_sel}.control_valsl)';
        
    end
           
    % derivatives wrt tf (phase A)
    
    if problem.structure{phase_sel}.imposed_tf == 0
        
        temp_tf = phase_boundaries{phase_sel}.tf;                           % get t0 variables of this phase
        temp_tf = temp_tf+h;                                                % add small increment
        
        phase_boundaries_temp = phase_boundaries;                           % clone phase_boundaries structure
        phase_boundaries_temp{phase_sel}.tf = temp_tf;                      % modify temporary clone, adding the incremented static variables
        
        c2 = impose_linkage(phase_boundaries_temp,this_phase_id,next_phase_id,problem,0);   % call this function, recursively, to get the values
        Jtf = (c2-c)/h;                                                     % compute finite difference approx of Jacobian
        
        Jc.tf_a = Jtf';
        
    end
    
    % derivatives wrt xf (phase A)
    
    ids = 1:length(phase_boundaries{phase_sel}.xf);                         % indexing vector for id of free xf
    ids = ids(~problem.structure{phase_sel}.imposed_final_states);          % keep only indexes of free xf
    nxf = length(ids);                                                      % number of free final states of this phase
    
    if nxf>0
        
        Jxf = zeros(size(c,1),nxf);
        
        for i=1:nxf
            
            temp_xf = phase_boundaries{phase_sel}.xf;                       % get final state variables of this phase
            temp_xf(ids(i)) = temp_xf(ids(i))+h;                            % add small increment
            
            phase_boundaries_temp = phase_boundaries;                       % clone phase_boundaries structure
            phase_boundaries_temp{phase_sel}.xf = temp_xf;                  % modify temporary clone, adding the incremented static variables
            
            c2 = impose_linkage(phase_boundaries_temp,this_phase_id,next_phase_id,problem,0);   % call this function, recursively, to get the values
            Jxf(:,i) = (c2-c)'/h;                                            % compute finite difference approx of Jacobian
            
        end
        
        Jc.xf_a = Jxf';
        
    end
    
    % derivatives wrt uf (phase A)
    
    nuf = length(phase_boundaries{phase_sel}.uf);                           % get number of control variables of this phase
    
    if nuf>0
        
        Juf = zeros(size(c,1),nuf);
        
        for i=1:nuf
            
            temp_uf = phase_boundaries{phase_sel}.uf;                       % get final control VALUES of this phase
            temp_uf(i) = temp_uf(i)+h;                                      % add small increment
            
            phase_boundaries_temp = phase_boundaries;                       % clone phase_boundaries structure
            phase_boundaries_temp{phase_sel}.uf = temp_uf;                  % modify temporary clone, adding the incremented static variables
            
            c2 = impose_linkage(phase_boundaries_temp,this_phase_id,next_phase_id,problem,0);   % call this function, recursively, to get the values
            Juf(:,i) = (c2-c)'/h;                                        % compute finite difference approx of Jacobian
            
        end
        
        % SUPER IMPORTANT: NEED TO CAST THIS DERIVATIVE IN THE SPACE OF
        % THE ACTUAL VARIABLES, SO WE NEED TO MULTIPLY THOSE
        % DERIVATIVES BY THE BASIS FUNCTIONS OF THE CONTROLS (LIKE FOR
        % THE DERIVATIVES OF THE OTHER CONSTRAINTS)
        
        Jc.uf_a = (Juf*problem.structure{phase_sel}.control_valsr)';
        
    end
    
    %% Derivatives wrt phase B variables
    
    phase_sel = next_phase_id;

    % derivatives wrt static variables (phase B)
    
    nstatic_b = length(phase_boundaries{phase_sel}.static);                 % get number of static variables of this phase
    
    if nstatic_b>0
        
        Jstatic = zeros(size(c,1),nstatic_b);
        
        for i=1:nstatic_b
            
            temp_static = phase_boundaries{phase_sel}.static;               % get static variables of this phase
            temp_static(i) = temp_static(i)+h;                              % add small increment
            
            phase_boundaries_temp = phase_boundaries;                       % clone phase_boundaries structure
            phase_boundaries_temp{phase_sel}.static = temp_static;          % modify temporary clone, adding the incremented static variables
            
            c2 = impose_linkage(phase_boundaries_temp,this_phase_id,next_phase_id,problem,0);   % call this function, recursively, to get the values
            Jstatic(:,i) = (c2-c)'/h;                                       % compute finite difference approx of Jacobian
            
        end
        
        Jc.static_b = Jstatic';                                             % fmincon uses the Transposed of the Jacobian...
        
    end
    
    % derivatives wrt t0 (phase B)
    
    if problem.structure{phase_sel}.imposed_t0 == 0
        
        temp_t0 = phase_boundaries{phase_sel}.t0;                           % get t0 variables of this phase
        temp_t0 = temp_t0+h;                                                % add small increment
        
        phase_boundaries_temp = phase_boundaries;                           % clone phase_boundaries structure
        phase_boundaries_temp{phase_sel}.t0 = temp_t0;                      % modify temporary clone, adding the incremented static variables
        
        c2 = impose_linkage(phase_boundaries_temp,this_phase_id,next_phase_id,problem,0);   % call this function, recursively, to get the values
        Jt0 = (c2-c)/h;                                                     % compute finite difference approx of Jacobian
        
        Jc.t0_b = Jt0';
        
    end
    
    % derivatives wrt x0 (phase B)
    
    ids = 1:length(phase_boundaries{phase_sel}.x0);                         % indexing vector for id of free x0
    ids = ids(~problem.structure{phase_sel}.imposed_initial_states);        % keep only indexes of free x0
    nx0 = length(ids);                                                      % number of free initial states of this phase
    
    if nx0>0
        
        Jx0 = zeros(size(c,1),nx0);
        
        for i=1:nx0
            
            temp_x0 = phase_boundaries{phase_sel}.x0;                       % get initial state variables of this phase
            temp_x0(ids(i)) = temp_x0(ids(i))+h;                            % add small increment
            
            phase_boundaries_temp = phase_boundaries;                       % clone phase_boundaries structure
            phase_boundaries_temp{phase_sel}.x0 = temp_x0;                  % modify temporary clone, adding the incremented static variables
            
            c2 = impose_linkage(phase_boundaries_temp,this_phase_id,next_phase_id,problem,0);   % call this function, recursively, to get the values
            Jx0(:,i) = (c2-c)'/h;                                            % compute finite difference approx of Jacobian
            
        end
        
        Jc.x0_b = Jx0';
        
    end
    
    % derivatives wrt u0 (phase B)
    
    nu0 = length(phase_boundaries{phase_sel}.u0);                           % get number of control variables in this phase
    
    if nu0>0
        
        Ju0 = zeros(size(c,1),nu0);
        
        for i=1:nu0
            
            temp_u0 = phase_boundaries{phase_sel}.u0;                       % get final control VALUES of this phase
            temp_u0(i) = temp_u0(i)+h;                                      % add small increment
            
            phase_boundaries_temp = phase_boundaries;                       % clone phase_boundaries structure
            phase_boundaries_temp{phase_sel}.u0 = temp_u0;                  % modify temporary clone, adding the incremented static variables
            
            c2 = impose_linkage(phase_boundaries_temp,this_phase_id,next_phase_id,problem,0);   % call this function, recursively, to get the values
            Ju0(:,i) = (c2-c)'/h;                                        % compute finite difference approx of Jacobian
            
        end
        
        % SUPER IMPORTANT: NEED TO CAST THIS DERIVATIVE IN THE SPACE OF
        % THE ACTUAL VARIABLES, SO WE NEED TO MULTIPLY THOSE
        % DERIVATIVES BY THE BASIS FUNCTIONS OF THE CONTROLS (LIKE FOR
        % THE DERIVATIVES OF THE OTHER CONSTRAINTS)
        
        Jc.u0_b = (Ju0*problem.structure{phase_sel}.control_valsl)';
        
    end
           
    % derivatives wrt tf (phase B)
    
    if problem.structure{phase_sel}.imposed_tf == 0
        
        temp_tf = phase_boundaries{phase_sel}.tf;                           % get t0 variables of this phase
        temp_tf = temp_tf+h;                                                % add small increment
        
        phase_boundaries_temp = phase_boundaries;                           % clone phase_boundaries structure
        phase_boundaries_temp{phase_sel}.tf = temp_tf;                      % modify temporary clone, adding the incremented static variables
        
        c2 = impose_linkage(phase_boundaries_temp,this_phase_id,next_phase_id,problem,0);   % call this function, recursively, to get the values
        Jtf = (c2-c)/h;                                                     % compute finite difference approx of Jacobian
        
        Jc.tf_b = Jtf';
        
    end
    
    % derivatives wrt xf (phase B)
    
    ids = 1:length(phase_boundaries{phase_sel}.xf);                         % indexing vector for id of free xf
    ids = ids(~problem.structure{phase_sel}.imposed_final_states);          % keep only indexes of free xf
    nxf = length(ids);                                                      % number of free final states of this phase
    
    if nxf>0
        
        Jxf = zeros(size(c,1),nxf);
        
        for i=1:nxf
            
            temp_xf = phase_boundaries{phase_sel}.xf;                       % get final state variables of this phase
            temp_xf(ids(i)) = temp_xf(ids(i))+h;                            % add small increment
            
            phase_boundaries_temp = phase_boundaries;                       % clone phase_boundaries structure
            phase_boundaries_temp{phase_sel}.xf = temp_xf;                  % modify temporary clone, adding the incremented static variables
            
            c2 = impose_linkage(phase_boundaries_temp,this_phase_id,next_phase_id,problem,0);   % call this function, recursively, to get the values
            Jxf(:,i) = (c2-c)'/h;                                            % compute finite difference approx of Jacobian
            
        end
        
        Jc.xf_b = Jxf';
        
    end
    
    % derivatives wrt uf (phase B)
    
    nuf = length(phase_boundaries{phase_sel}.uf);                           % get number of control variables of this phase
    
    if nuf>0
        
        Juf = zeros(size(c,1),nuf);
        
        for i=1:nuf
            
            temp_uf = phase_boundaries{phase_sel}.uf;                       % get final control VALUES of this phase
            temp_uf(i) = temp_uf(i)+h;                                      % add small increment
            
            phase_boundaries_temp = phase_boundaries;                       % clone phase_boundaries structure
            phase_boundaries_temp{phase_sel}.uf = temp_uf;                  % modify temporary clone, adding the incremented static variables
            
            c2 = impose_linkage(phase_boundaries_temp,this_phase_id,next_phase_id,problem,0);   % call this function, recursively, to get the values
            Juf(:,i) = (c2-c)'/h;                                        % compute finite difference approx of Jacobian
            
        end
        
        % SUPER IMPORTANT: NEED TO CAST THIS DERIVATIVE IN THE SPACE OF
        % THE ACTUAL VARIABLES, SO WE NEED TO MULTIPLY THOSE
        % DERIVATIVES BY THE BASIS FUNCTIONS OF THE CONTROLS (LIKE FOR
        % THE DERIVATIVES OF THE OTHER CONSTRAINTS)
        
        Jc.uf_b = (Juf*problem.structure{phase_sel}.control_valsr)';
        
    end
    
    %% Cast Jacobians as sparse variables to help fmincon

    Jc.static_a = sparse(Jc.static_a);    
    Jc.t0_a = sparse(Jc.t0_a);
    Jc.x0_a = sparse(Jc.x0_a);
    Jc.u0_a = sparse(Jc.u0_a);        
    Jc.tf_a = sparse(Jc.tf_a);
    Jc.xf_a = sparse(Jc.xf_a);
    Jc.uf_a = sparse(Jc.uf_a);    
 
    Jc.static_b = sparse(Jc.static_b);    
    Jc.t0_b = sparse(Jc.t0_b);
    Jc.x0_b = sparse(Jc.x0_b);
    Jc.u0_b = sparse(Jc.u0_b);        
    Jc.tf_b = sparse(Jc.tf_b);
    Jc.xf_b = sparse(Jc.xf_b);
    Jc.uf_b = sparse(Jc.uf_b);    
    
end

end
