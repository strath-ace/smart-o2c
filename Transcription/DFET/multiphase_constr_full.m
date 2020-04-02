function [c_all,ceq_all,Jc_all,Jceq_all] = multiphase_constr_full (problem,x_in,jacflag,tstart,tlim)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%
% This function computes all constraints and interphase linkage equations
% and their Jacobians, for all phases.

%% Extract state, control and static variables of all phases

x_in = x_in(:); % ensure column vector, bad bad things happen without warning, otherwise...

c_all = [];
ceq_all = [];
Jc_all = [];
Jceq_all = [];

x_guess_phase = cell(problem.num_phases,1);

phase_boundaries = cell(problem.num_phases,1);

for i = 1:problem.num_phases
    
    x_guess_phase{i} = x_in(problem.phase_mask==i);
    
    % Get x0, xf, u0, uf, t0, tf and static variables of this phase
    
    [x_0,u_0,t_0,x_f,u_f,t_f,static] = get_phase_boundary_vals (x_guess_phase{i},problem.structure{i});
    
    phase_boundaries{i}.x0 = x_0;
    phase_boundaries{i}.u0 = u_0;
    phase_boundaries{i}.t0 = t_0;
    phase_boundaries{i}.xf = x_f;
    phase_boundaries{i}.uf = u_f;
    phase_boundaries{i}.tf = t_f;
    phase_boundaries{i}.static = static;
    phase_boundaries{i}.constants = problem.structure{i}.constants;
    phase_boundaries{i}.scales = problem.structure{i}.scales;
    
end

%% Compute all constraints and Jacobians of all phases

Jc = cell(problem.num_phases,1);
Jceq = cell(problem.num_phases,1);
Jclink = cell(problem.num_phases);
clink = cell(problem.num_phases);

for i = 1:problem.num_phases
    
    % Compute all constraints of phase i
    [c,ceq,Jc{i},Jceq{i}] = Full_constr_dynamics (problem.structure{i},x_guess_phase{i},jacflag);
    
    % Append constraints of this phase to existing constraints
    
    c_all = [c_all;c];
    ceq_all = [ceq_all;ceq];
    
    % Compute linkage constraints (if next phase is defined for this phase)
    
    if ~isempty(problem.structure{i}.next_phases)
        
        nlinks = length(problem.structure{i}.next_phases);

        for j = 1:nlinks
            
            next_phase_id = problem.structure{i}.next_phases(j);
        
            [clink{i,next_phase_id},Jclink{i,next_phase_id}] = impose_linkage(phase_boundaries,i,next_phase_id,problem,jacflag); % Jclink{i}{j} contains Jacobian of linkage function between phase i and phase j
        
            ceq_all = [ceq_all;clink{i,next_phase_id}];
            
        end
        
    end
    
end

%% Assemble Jacobians

if jacflag == 1
    
    startx_eq = 1;
    startx_ineq = 1;
    %starty_eq = 1;
    %starty_ineq = 1;
    
    for i = 1:problem.num_phases

        % Jacobians of equality and inequality constraints of phase i
        
        endx_eq = startx_eq+size(Jceq{i},2)-1;
        endx_ineq = startx_ineq+size(Jc{i},2)-1;
        %endy_eq = starty_eq+size(Jceq{i},1)-1;
        %endy_ineq = starty_ineq+size(Jc{i},1)-1;

        position = logical((problem.phase_mask==i));
        Jceq_all(position,startx_eq:endx_eq) = Jceq{i};
        Jc_all(position,startx_ineq:endx_ineq) = Jc{i};
    
        startx_eq = endx_eq+1;
        startx_ineq = endx_ineq+1;
            
        if ~isempty(problem.structure{i}.next_phases)
       
            for links = 1:length(problem.structure{i}.next_phases)
                
                j = problem.structure{i}.next_phases(links);          
                
                endx_eq = startx_eq+size(clink{i,j},1)-1;
                
                %% Jacobians wrt vars of phase i
            
                % Jacobians of linkage function between phase i and j wrt
                % static vars of phase i (static_a)

                position = logical((problem.phase_mask==i).*problem.other_vars);
                Jceq_all(position,startx_eq:endx_eq) = Jclink{i,j}.static_a;

                % Jacobians of linkage function between phase i and j wrt
                % u0 of phase i (u0_a)

                position = logical((problem.phase_mask==i).*problem.u0_vars);
                Jceq_all(position,startx_eq:endx_eq) = Jclink{i,j}.u0_a;

                % Jacobians of linkage function between phase i and j wrt
                % x0 vars of phase i (x0_a)

                position = logical((problem.phase_mask==i).*problem.x0_vars);
                Jceq_all(position,startx_eq:endx_eq) = Jclink{i,j}.x0_a;

                % Jacobians of linkage function between phase i and j wrt
                % t0 vars of phase i (t0_a)

                position = logical((problem.phase_mask==i).*problem.t0_vars);
                Jceq_all(position,startx_eq:endx_eq) = Jclink{i,j}.t0_a;
                
                % Jacobians of linkage function between phase i and j wrt
                % uf vars of phase i (uf_a)

                position = logical((problem.phase_mask==i).*problem.uf_vars);
                Jceq_all(position,startx_eq:endx_eq) = Jclink{i,j}.uf_a;

                % Jacobians of linkage function between phase i and j wrt
                % xf vars of phase i (xf_a)

                position = logical((problem.phase_mask==i).*problem.xf_vars);
                Jceq_all(position,startx_eq:endx_eq) = Jclink{i,j}.xf_a;

                % Jacobians of linkage function between phase i and j wrt
                % tf of phase i (tf_a)

                position = logical((problem.phase_mask==i).*problem.tf_vars);
                Jceq_all(position,startx_eq:endx_eq) = Jclink{i,j}.tf_a;

                %% Jacobians wrt vars of phase j
            
                % Jacobians of linkage function between phase i and j wrt
                % static vars of phase j (static_b)

                position = logical((problem.phase_mask==j).*problem.other_vars);
                Jceq_all(position,startx_eq:endx_eq) = Jclink{i,j}.static_b;

                % Jacobians of linkage function between phase i and j wrt
                % u0 vars of phase j (u0_b)

                position = logical((problem.phase_mask==j).*problem.u0_vars);
                Jceq_all(position,startx_eq:endx_eq) = Jclink{i,j}.u0_b;

                % Jacobians of linkage function between phase i and j wrt
                % x0 vars of phase j (x0_b)

                position = logical((problem.phase_mask==j).*problem.x0_vars);
                Jceq_all(position,startx_eq:endx_eq) = Jclink{i,j}.x0_b;

                % Jacobians of linkage function between phase i and j wrt
                % t0 vars of phase j (t0_b)

                position = logical((problem.phase_mask==j).*problem.t0_vars);
                Jceq_all(position,startx_eq:endx_eq) = Jclink{i,j}.t0_b;
                
                % Jacobians of linkage function between phase i and j wrt
                % uf vars of phase j (uf_b)

                position = logical((problem.phase_mask==j).*problem.uf_vars);
                Jceq_all(position,startx_eq:endx_eq) = Jclink{i,j}.uf_b;

                % Jacobians of linkage function between phase i and j wrt
                % xf vars of phase j (xf_b)

                position = logical((problem.phase_mask==j).*problem.xf_vars);
                Jceq_all(position,startx_eq:endx_eq) = Jclink{i,j}.xf_b;

                % Jacobians of linkage function between phase i and j wrt
                % tf vars of phase j (tf_b)

                position = logical((problem.phase_mask==j).*problem.tf_vars);
                Jceq_all(position,startx_eq:endx_eq) = Jclink{i,j}.tf_b;
                
                startx_eq = endx_eq+1;

            end
            
        end
           
    end
    
    % cast all into a sparse matrix to simplify the work of fmincon
    
    Jceq_all = sparse(Jceq_all);
    Jc_all = sparse(Jc_all);
    
end

%% Check for NANs or INFs

% if isnan(Jceq_all)
%     
%     keyboard
%     
% end
% 
% if isnan(Jc_all)
%     
%     keyboard
%     
% end
% 
% if isinf(Jceq_all)
%     
%     keyboard
%     
% end
% 
% if isinf(Jc_all)
%     
%     keyboard
%     
% end
% 
% if isnan(ceq_all)
%     
%     keyboard
%     
% end
% 
% if isnan(c_all)
%     
%     keyboard
%     
% end
% 
% if isinf(ceq_all)
%     
%     keyboard
%     
% end
% 
% if isinf(c_all)
%     
%     keyboard
%     
% end

%% Check for timeout

if tlim>0  % allows to give infinite time by setting a negative limit

    tnow = toc(tstart);

    if tnow>tlim % timeout
    
        errorStruct.message = ['Timeout of ',num2str(tlim),' seconds reached'];
        errorStruct.identifier = 'multiphase_constr_full:timeOutReached';
        error(errorStruct);

    end    
    
end


end