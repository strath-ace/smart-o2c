function [J_out,C,Ceq,GJ,GC,GCeq,ToF,Equin_sorted,Equin_sorted_ct,time_sorted,J,E0] = ...
    FABLE_transcription_MS(x, parameters, options, constants)

% =========================================================================
% Multiple shooting transcription of the low-thrust problem for FABLE.
%
% Function for the computation of the objective function
% and of the non-linear equality constraints.
% =========================================================================
% Input: x          -> solution vector
%        parameters -> structure containing the parameters of the problem
%        options    -> structure containing the options for the problem
%        constants  -> structure of constants

% Output: J_out     -> objective function
%         C         -> non-linear inequality constraints
%         Ceq       -> non-linear equalities constraints
%         GCeq      -> gradients of the Ceq
%         GJ        -> gradients of J
%         ToF       -> time of flight for the transfer
%         Equin_sorted -> cells with solution data. Each cell contains the
%                         equinoctial elements and mass variation of each
%                         thrust or coast arc
%         Equin_sorted_ct -> vectors of 0s and 1s (1 for thurst arcs, 0 for
%                            coast arcs)
%         time_sorted -> same as Equin_sorted, but for time variation
%         J           ->
% For FABLE1:
% The vector x has 1x14*arcs variables, where arcs is the number of thrust arcs.
% For each arc, 14 variable are defined:
% - alpha (azimuth)
% - beta(elevation)
% - 6 equinoctial variables for the initial point of the propulsed arc
% - 6 equinoctial variabes for the final point of the propulsed arc.
% They are organized in the x vector as follows:
% [alpha, beta, [a, P1, P2, Q1, Q2, L], [a, P1, P2, Q1, Q2, Delta L]]
% for each propulsed arc.

% For FABLE2:
% The vector x has 8*arcs+3 variables, where arcs is the number of arcs.
% For each (non-initial) arc, 8 variables are defined:
% - epsilon (acceleration)
% - alpha (azimuth)
% - beta(elevation)
% - 5 equinoctial variables for the initial point of the arc
% 3 variables are used for the first arc (epsilon, azimuth and elevation -
% 5 equinoctial elements are given in this case by the departure
% conditions)
% They are organized in the x vector as follows:
% [[epsilon, alpha, beta],[epsilon, alpha, beta], [a, P1, P2, Q1, Q2], ...
% [epsilon, alpha, beta], [a, P1, P2, Q1, Q2], [epsilon, alpha, beta], [a,
% P1, P2, Q1, Q2]....]
% for each arc.
% =========================================================================
% Author: Marilena Di Carlo
% email: marilena.di-carlo@strath.ac.uk
% =========================================================================

%% Initialisation

% Number of transfer arcs
arcs         = parameters.arcs;

% Equinoctial elements arrival point
arrival_eq   = parameters.arrival_eq;

% Time of flight
curr_tof     = parameters.curr_tof;

% Number of steps for analytic integration
n            = parameters.n;


% Equinoctial elements departure point
if options.FABLE1 && numel(x) == 14*arcs || ...
        options.FABLE2 && numel(x) == 8*(arcs-1)+3
    % If departure state is given and declination and azimuth at departure
    % do not have to be optimised
    departure_eq = parameters.departure_eq;
elseif options.FABLE1 && numel(x) == 14*arcs + 2 || ...
        options.FABLE2 && numel(x) == 8*(arcs-1)+3 + 2
    % Departure position
    departure_pos = parameters.dep_body_r;
    % Velocity at departure is given by the sum of the velocity of the
    % central body and the velocity vinf
    v_inf = parameters.vinf0/ constants.DU * constants.TU * constants.sec_day;
    % Azimuth and declination
    alpha = x(end-1);
    delta = x(end);
    
    departure_vel = parameters.dep_body_v + ...
        v_inf * [cos(alpha) * cos(delta), sin(alpha) * cos(delta), sin(delta)];
    
    departure_kep = cart2kep([departure_pos, departure_vel], constants.mu);
    
    departure_eq = kep2eq(departure_kep);
end

% Initial conditions: six equinoctial elements, mass
E0 = [departure_eq; parameters.m];


% Cost function
J = 0;

% Equality constraints
Ceq = [];

% Time of flight
ToF = 0;

% Initialize variables for coast and thrust states
Equin_coast = [];
Equin_thrust = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FABLE1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if options.FABLE1
    
    % Low-thrust acceleration
    eps_max      = parameters.eps_max;
    
    % Variable to define position of variables in vector xopt;
    k = numel(x) - 14*arcs;
    
    %% FIRST PROPULSED ARC
    
    % Initial equinoctial element for forward propagation on the first low
    % thrust arc
    Equin_initial_prop_forw  = x(3:8);
    
    % Final equinoctial elements of first forward propagation on the first low
    % thrust arc
    Equin_final_prop_forw = [x(9:13) x(8)+x(14)];
    
    % True longitude variation over first low thrust arc
    L_prop_forw  = linspace(Equin_initial_prop_forw(6),  Equin_final_prop_forw(6), n);
    
    
    %% LAST PROPULSED ARC
    
    % Initial equinoctial element for backward propagation on the first LT
    % backward propagated arc (last LT arcs)
    Equin_initial_prop_back  = [x(end-5-k:end-1-k) x(end-6-k)+x(end-k)];
    
    % Final equinoctial elements for backward propagation on the first LT
    % backward propagated arc (last LT arcs)
    Equin_final_prop_back = x(end-11-k:end-6-k);
    
    % Longitude values for propagation
    L_prop_back  = linspace(Equin_initial_prop_back(6),  Equin_final_prop_back(6), n);
    
    %% Initialise GCeq
    
    
    if options.const_ToF
        GCeq = zeros(5 * (2*arcs+1)  + 1, numel(x));
    else
        GCeq = zeros(5 * (2*arcs+1), numel(x));
    end
    
    t_arcs = [];
    t_arcs_var = [];
    t_arcs_ct = [];
    index_t = 1;
    
    delta_FD = 1e-8;
    
    %% FIRST COAST ARC?
    
    % Define an initial coast arc only if the first point of the vector x
    % does not coincide with the departure point
    
    if departure_eq(6) == Equin_initial_prop_forw(6)
        
        % If the longitude of the initial point of the first LT arc coincides with the
        % departure longitude, there is no need of using a coast arc and the
        % first equality contraint (coincidence of the two points) has to be
        % defined
        
        Ceq = [Ceq; departure_eq(1:5) - Equin_initial_prop_forw(1:5)'];
        
        
    else
        % If the initial point of the first propulsed arc does not coincide
        % with the departure point, use a coast arc to reach the first point
        
        % True longitude variation
        L_coast_forw = linspace(departure_eq(6), Equin_initial_prop_forw(6), n);
        
        % Propagation. Since this is a coast arc it is not necessary to choose
        % the model for the thrust. Therefore the following is commented and
        % substituted by
        [Equin_coast_forw, t] = AnEquin_forward_m(L_coast_forw, departure_eq, ...
            0, 0, 0, constants.mu);
        
        
        Equin_coast = [Equin_coast, Equin_coast_forw];
        
        Equin_sorted{1} = Equin_coast_forw;
        Equin_sorted_ct(1) = 0;
        
        time_sorted{1} = t;
        
        % Time of flight update
        ToF = ToF + t(end);
        t_arcs = [t_arcs; t(end)*ones(1,numel(x))];
        t_arcs_var = [t_arcs_var; t(end)*ones(1,numel(x))];
        t_arcs_ct = [t_arcs_ct; 0];
        
        index_t = index_t + 1;
        
        % The final point of the initial coast arc has to coincide with the initial point
        % of the first low thrust arc
        Ceq = [Ceq; Equin_coast_forw(1:5,end) - Equin_initial_prop_forw(1:5)'];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Variation of the true longitude of the 1st ON node
        if options.gradient == 1
            
            % Change of the true longitude of the 1st ON node
            x_var = x;
            x_var(8) = x_var(8) + delta_FD;
            
            % Coast arc
            eps_GCeq = 0;
            alpha_GCeq = 0;
            beta_GCeq = 0;
            % Final point propagation
            L_prop_GCeq = linspace(departure_eq(6), x_var(8), n);
            % Initial point propagation
            start_eq_GCeq = departure_eq;
            % Final state of propagation should be...:
            check_eq_GCeq = Equin_initial_prop_forw(1:5)';
            % Constraint under x_var
            [Ceq_var, t_var] = FABLE_grad_constr(eps_GCeq, alpha_GCeq, beta_GCeq, L_prop_GCeq, ...
                start_eq_GCeq, check_eq_GCeq, options, constants);
            % GCeq
            GCeq(1:5,8) = (Ceq_var - Ceq(1:5)) / delta_FD;
            
            %
            t_arcs_var(index_t - 1, 8) = t_var(end);
            
            % Variation of the azimuth and declination angles at departure
            if numel(x) == 14*arcs + 2
                
                % Change alpha and get effect on GCeq
                % Azimuth and declination
                alpha_var = x(end-1) + delta_FD;
                delta_var = x(end) + delta_FD;
                
                % ------------alpha variation
                departure_vel_var = parameters.dep_body_v + ...
                    v_inf * [cos(alpha_var) * cos(delta), sin(alpha_var) * cos(delta), sin(delta)];
                
                departure_kep_var = cart2kep([departure_pos, departure_vel_var], constants.mu);
                
                departure_eq_var = kep2eq(departure_kep_var);
                
                
                % Coast arc
                eps_GCeq = 0;
                alpha_GCeq = 0;
                beta_GCeq = 0;
                % Final point propagation
                L_prop_GCeq = linspace(departure_eq_var(6), Equin_initial_prop_forw(6), n);
                % Initial point propagation
                start_eq_GCeq = departure_eq_var;
                % Final state of propagation should be...:
                check_eq_GCeq = Equin_initial_prop_forw(1:5)';
                % Constraint under x_var
                [Ceq_var, t_var] = FABLE_grad_constr(eps_GCeq, alpha_GCeq, beta_GCeq, L_prop_GCeq, ...
                    start_eq_GCeq, check_eq_GCeq, options, constants);
                % GCeq
                GCeq(1:5,end-1) = (Ceq_var - Ceq(1:5)) / delta_FD;
                
                %
                t_arcs_var(index_t - 1, end-1) = t_var(end);
                
                
                
                
                % -----------delta variation
                departure_vel_var = parameters.dep_body_v + ...
                    v_inf * [cos(alpha) * cos(delta_var), sin(alpha) * cos(delta_var), sin(delta_var)];
                
                departure_kep_var = cart2kep([departure_pos, departure_vel_var], constants.mu);
                
                departure_eq_var = kep2eq(departure_kep_var);
                
                % Coast arc
                eps_GCeq = 0;
                alpha_GCeq = 0;
                beta_GCeq = 0;
                % Final point propagation
                L_prop_GCeq = linspace(departure_eq_var(6), Equin_initial_prop_forw(6), n);
                % Initial point propagation
                start_eq_GCeq = departure_eq_var;
                % Final state of propagation should be...:
                check_eq_GCeq = Equin_initial_prop_forw(1:5)';
                % Constraint under x_var
                [Ceq_var, t_var] = FABLE_grad_constr(eps_GCeq, alpha_GCeq, beta_GCeq, L_prop_GCeq, ...
                    start_eq_GCeq, check_eq_GCeq, options, constants);
                % GCeq
                GCeq(1:5,end) = (Ceq_var - Ceq(1:5)) / delta_FD;
                
                %
                t_arcs_var(index_t - 1, end) = t_var(end);
            end
            
            
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if options.gradient
        
        GCeq(1,3) = -1;
        GCeq(2,4) = -1;
        GCeq(3,5) = -1;
        GCeq(4,6) = -1;
        GCeq(5,7) = -1;
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% LAST COAST ARC?
    
    
    if arrival_eq(6) == Equin_initial_prop_back(6)
        
        Ceq = [Ceq; arrival_eq(1:5) - Equin_initial_prop_back(1:5)'];
        
    else
        
        L_coast_back = linspace(arrival_eq(6), Equin_initial_prop_back(6), n);
        
        [Equin_coast_back, t] = AnEquin_forward_m(L_coast_back, arrival_eq, ...
            0, 0, 0, constants.mu);
        
        
        Equin_coast = [Equin_coast, Equin_coast_back];
        
        
        Equin_sorted{2*arcs+1} = Equin_coast_back(:,end:-1:1);
        Equin_sorted_ct(2*arcs+1) = 0;
        
        time_sorted{2*arcs+1} = -t;
        
        
        
        % Time of flight update (negative time for backward propagation)
        ToF = ToF - t(end);
        t_arcs = [t_arcs; -t(end)*ones(1,numel(x))];
        t_arcs_var = [t_arcs_var; -t(end)*ones(1,numel(x))];
        t_arcs_ct = [t_arcs_ct; 0];
        index_t = index_t + 1;
        
        % Final point of the initial coast arc to coincide with the initial point
        % of the first backward propagated LT arc
        Ceq = [Ceq; Equin_coast_back(1:5,end) - Equin_initial_prop_back(1:5)'];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Variation of the true longitude of the last OFF node
        if options.gradient == 1
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % NEW TO COMMENT
            % Change LE
            index_GCeq = length(x) - 6-k;
            x_var = x;
            
            x_var(index_GCeq) = x_var(index_GCeq) + delta_FD;
            
            % Coast arc conditions
            eps_GCeq = 0;
            alpha_GCeq = 0;
            beta_GCeq = 0;
            % True longitudes for propagation
            L_prop_GCeq = linspace(arrival_eq(6), x_var(end-6-k)+x_var(end-k), n);
            % Initial point of propagation
            start_eq_GCeq = arrival_eq;
            % Desired final state of propagation
            check_eq_GCeq = Equin_initial_prop_back(1:5)';
            % Propagate single arc and compute variation in Ceq
            [~, t_var] = FABLE_grad_constr(eps_GCeq,  ...
                alpha_GCeq, beta_GCeq, L_prop_GCeq, start_eq_GCeq, check_eq_GCeq, options, constants);
            
            %
            t_arcs_var(index_t - 1, end-6-k) = -t_var(end);
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Gradient constraints
            x_var = x;
            x_var(end-k) = x_var(end-k) + delta_FD;
            
            % Coast arc conditions
            eps_GCeq = 0;
            alpha_GCeq = 0;
            beta_GCeq = 0;
            % True longitudes for propagation
            L_prop_GCeq = linspace(arrival_eq(6), x_var(end-6-k)+x_var(end-k), n);
            % Initial point of propagation
            start_eq_GCeq = arrival_eq;
            % Desired final state of propagation
            check_eq_GCeq = Equin_initial_prop_back(1:5)';
            % Propagate single arc and compute variation in Ceq
            [Ceq_var, t_var] = FABLE_grad_constr(eps_GCeq,  ...
                alpha_GCeq, beta_GCeq, L_prop_GCeq, start_eq_GCeq, check_eq_GCeq, options, constants);
            
            GCeq(6:10,end-k) = (Ceq_var - Ceq(6:10)) / delta_FD;
            
            %
            t_arcs_var(index_t - 1, end-k) = -t_var(end);
            
            
            
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end
    
    if options.gradient
        
        GCeq(6,end-k-5) = -1;
        GCeq(7,end-k-4) = -1;
        GCeq(8,end-k-3) = -1;
        GCeq(9,end-k-2) = -1;
        GCeq(10,end-k-1) = -1;
        
    end
    
    temp = 11;
    
    
    %%
    % Pre-allocation
    J_var = zeros(arcs,numel(x));
    J_all = [];
    
    % Equin_thrust = zeros(6, arcs * parameters.n);
    
    %%
    
    for i = 1 : arcs/2
        
        
        % ---------------------------------------------------------------------
        % Forward propagated low-thrust arc
        % ---------------------------------------------------------------------
        if options.fun == 0
            [EquinLT_forward, t] = AnEquin_forward_m(L_prop_forw,...
                Equin_initial_prop_forw', eps_max, x(1+(i-1)*14), x(2+(i-1)*14), ...
                constants.mu);
        elseif options.fun == 1
            [EquinLT_forward, t] = AnEquin_forward_m_r2(L_prop_forw,...
                Equin_initial_prop_forw', eps_max, x(1+(i-1)*14), x(2+(i-1)*14), ...
                constants.mu);
        end
        
        
        Equin_thrust = [Equin_thrust, EquinLT_forward];
        
        Equin_sorted{2+(i-1)*2} = EquinLT_forward;
        Equin_sorted_ct(2+(i-1)*2) = 1;
        
        time_sorted{2+(i-1)*2} = t;
        
        if any(time_sorted{2+(i-1)*2})<0
            keyboard
        end
        
        if options.fun == 0
            % Cost function: time spent using propulsion
            J = J + t(end);
        elseif options.fun == 1
            % Cost function: analytic expression DeltaV
            B = sqrt(1 - Equin_initial_prop_forw(2)^2 - Equin_initial_prop_forw(3)^2);
            J = J + eps_max / (B * sqrt(constants.mu * Equin_initial_prop_forw(1)) )* ...
                (L_prop_forw(end) - L_prop_forw(1));
            J_all = [J_all; eps_max / (B * sqrt(constants.mu * Equin_initial_prop_forw(1)) )* ...
                (L_prop_forw(end) - L_prop_forw(1))];
        end
        
        % Time of flight update
        ToF = ToF + t(end);
        t_arcs = [t_arcs; t(end)*ones(1,numel(x))];
        t_arcs_var = [t_arcs_var; t(end)*ones(1,numel(x))];
        t_arcs_ct = [t_arcs_ct; 1];
        index_t = index_t + 1;
        
        % Matching condition between final point of forward propagated propulsed
        % arc and required final point of forward propagated propulsed arc
        Ceq = [Ceq; EquinLT_forward(1:5,end) - Equin_final_prop_forw(1:5)'];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if options.gradient == 1
            % Index corresponding to the position in the vector x that can cause
            % variation to the considered constraints
            index_GCeq = [1 2 3 4 5 6 7 8 14] + (i-1)*14;
            index_GJ = [3 4 5 14] + (i-1)*14;
            
            % Acceleration
            eps_GCeq = eps_max;
            
            for i_GCeq = 1 : length(index_GCeq)
                
                % Change one single variable in the vector x by delta
                x_var = x;
                x_var(index_GCeq(i_GCeq)) = x_var(index_GCeq(i_GCeq)) + delta_FD;
                
                % Azimuth and elevation angles
                alpha_GCeq = x_var(1+(i-1)*14);
                beta_GCeq  = x_var(2+(i-1)*14);
                
                % Initial point propagation
                start_eq_GCeq = x_var(3+14*(i-1): 8+14*(i-1));
                
                % Final state of propgation should be...
                check_eq_GCeq = [x_var(9+14*(i-1): 13+14*(i-1))'; ...
                    x_var(8+14*(i-1))+x_var(14+14*(i-1))];
                
                % True longitude variations
                L_prop_GCeq = linspace(start_eq_GCeq(6), check_eq_GCeq(6), n);
                
                [Ceq_var, t_var] = FABLE_grad_constr(eps_GCeq, alpha_GCeq, beta_GCeq,...
                    L_prop_GCeq, start_eq_GCeq, check_eq_GCeq(1:5), options, constants);
                
                GCeq(temp : temp+4, index_GCeq(i_GCeq)) = ...
                    (Ceq_var - Ceq(temp:temp+4)) / delta_FD;
                
                t_arcs_var(index_t - 1, index_GCeq(i_GCeq)) = t_var(end);
                
                % GJ for acceleration as 1/r2
                if ismember(index_GCeq(i_GCeq), index_GJ) && options.fun == 1
                    B_var = sqrt(1 - x_var(4 + (i-1)*14)^2 - ...
                        x_var(5 + (i-1)*14)^2);
                    J_var(1+2*(i-1),index_GCeq(i_GCeq)) = J_var(1+2*(i-1),index_GCeq(i_GCeq)) + ...
                        eps_max / (B_var * sqrt(constants.mu * x_var(3 + (i-1)*14)) )* ...
                        (x_var(14 + (i-1)*14));
                end
            end
            
            GCeq(temp,   9 + (i-1)*14) = -1;
            GCeq(temp+1, 10 + (i-1)*14) = -1;
            GCeq(temp+2, 11 + (i-1)*14) = -1;
            GCeq(temp+3, 12 + (i-1)*14) = -1;
            GCeq(temp+4, 13 + (i-1)*14) = -1;
            
            temp = temp + 5;
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % ---------------------------------------------------------------------
        % Forward coast arc from final point of low thurst arc to initial point
        % of next propulsed arc
        % ---------------------------------------------------------------------
        
        % Forward coast arc
        L_coast_forw = linspace(Equin_final_prop_forw(6), x(22+14*(i-1)), n);
        Equin_initial_coast_forw = Equin_final_prop_forw';
        
        [Equin_coast_forward, t] = AnEquin_forward_m(L_coast_forw, ...
            Equin_initial_coast_forw, 0, 0, 0, constants.mu);
        
        Equin_coast = [Equin_coast, Equin_coast_forward];
        
        Equin_sorted{3+(i-1)*2} = Equin_coast_forward;
        Equin_sorted_ct(3+(i-1)*2) = 0;
        
        time_sorted{3+(i-1)*2} = t;
        
        if any(time_sorted{3+(i-1)*2})<0
            keyboard
        end
        
        % Time of flight update
        ToF = ToF + t(end);
        t_arcs = [t_arcs; t(end)*ones(1,numel(x))];
        t_arcs_var = [t_arcs_var; t(end)*ones(1,numel(x))];
        t_arcs_ct = [t_arcs_ct; 0];
        index_t = index_t + 1;
        
        % Matching condition between final point of forward coast arc and
        % initial point of next propulsed arc
        
        Ceq = [Ceq; Equin_coast_forward(1:5,end) - x(17+14*(i-1):21+14*(i-1))'];
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if options.gradient == 1
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % NEW TO COMMENT
            index_GCeq = 8+ (i-1)*14;
            x_var = x;
            
            x_var(index_GCeq) = x_var(index_GCeq) + delta_FD;
            eps_GCeq = 0;
            alpha_GCeq = 0;
            beta_GCeq  = 0;
            
            % Initial point propagation
            start_eq_GCeq = [x_var(9+14*(i-1): 13+14*(i-1)) x_var(8+14*(i-1))+x_var(14+14*(i-1))];
            
            % Final desired state
            check_eq_GCeq = x_var(17 + 14 * (i-1) : 22+14*(i-1))';
            
            % True longitude variations
            L_prop_GCeq = linspace(start_eq_GCeq(6), check_eq_GCeq(6), n);
            
            %
            [Ceq_var, t_var] = FABLE_grad_constr(eps_GCeq, alpha_GCeq, beta_GCeq,...
                L_prop_GCeq, start_eq_GCeq, check_eq_GCeq(1:5), options, constants);
            
            t_arcs_var(index_t - 1, index_GCeq) = t_var(end);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            % Gradient constraints
            index_GCeq =  [9 10 11 12 13 14 22] + (i-1)*14;
            
            
            
            % Coast arcs
            eps_GCeq = 0;
            
            for i_GCeq = 1 : length(index_GCeq)
                
                x_var = x;
                
                x_var(index_GCeq(i_GCeq)) = x_var(index_GCeq(i_GCeq)) + delta_FD;
                
                alpha_GCeq = 0;
                beta_GCeq  = 0;
                
                % Initial point propagation
                start_eq_GCeq = [x_var(9+14*(i-1): 13+14*(i-1)) x_var(8+14*(i-1))+x_var(14+14*(i-1))];
                
                % Final desired state
                check_eq_GCeq = x_var(17 + 14 * (i-1) : 22+14*(i-1))';
                
                % True longitude variations
                L_prop_GCeq = linspace(start_eq_GCeq(6), check_eq_GCeq(6), n);
                
                %
                [Ceq_var, t_var] = FABLE_grad_constr(eps_GCeq, alpha_GCeq, beta_GCeq,...
                    L_prop_GCeq, start_eq_GCeq, check_eq_GCeq(1:5), options, constants);
                
                GCeq(temp : temp + 4,index_GCeq(i_GCeq)) = ...
                    (Ceq_var - Ceq(temp:temp+4)) / delta_FD;
                
                t_arcs_var(index_t - 1, index_GCeq(i_GCeq)) = t_var(end);
                
            end
            
            GCeq(temp,   17 + (i-1)*14) = -1;
            GCeq(temp+1, 18 + (i-1)*14) = -1;
            GCeq(temp+2, 19 + (i-1)*14) = -1;
            GCeq(temp+3, 20 + (i-1)*14) = -1;
            GCeq(temp+4, 21 + (i-1)*14) = -1;
            
            temp = temp + 5;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ---------------------------------------------------------------------
        % Backward propagated low-thrust arc
        % ---------------------------------------------------------------------
        
        if options.fun == 0
            [EquinLT_backward, t] = AnEquin_forward_m(L_prop_back, ...
                Equin_initial_prop_back', eps_max, x(end-13-k-14*(i-1)), ...
                x(end-12-k-14*(i-1)), constants.mu);
        elseif options.fun == 1
            [EquinLT_backward, t] = AnEquin_forward_m_r2(L_prop_back, ...
                Equin_initial_prop_back', eps_max, x(end-13-k-14*(i-1)), ...
                x(end-12-k-14*(i-1)), constants.mu);
        end
        
        Equin_thrust = [Equin_thrust, EquinLT_backward];
        
        
        Equin_sorted{2*arcs-(i-1)*2} = EquinLT_backward(:,end:-1:1);
        Equin_sorted_ct(2*arcs-(i-1)*2) = 1;
        
        time_sorted{2*arcs-(i-1)*2} = -t;
        if any(time_sorted{2*arcs-(i-1)*2})<0
            keyboard
        end
        
        if options.fun ==0
            % Cost function: time spent using propulsion
            J = J - t(end);
        elseif options.fun == 1
            % Cost function: analytic computation DeltaV
            B = sqrt(1 - Equin_initial_prop_back(2)^2 - Equin_initial_prop_back(3)^2);
            J = J + eps_max / (B * sqrt(constants.mu * Equin_initial_prop_back(1)) )* ...
                (L_prop_back(1) - L_prop_back(end));
            J_all = [J_all; eps_max / (B * sqrt(constants.mu * Equin_initial_prop_back(1)) )* ...
                (L_prop_back(1) - L_prop_back(end))];
        end
        
        
        % Time of flight update
        ToF = ToF - t(end);
        t_arcs = [t_arcs; -t(end)*ones(1,numel(x))];
        t_arcs_var = [t_arcs_var; -t(end)*ones(1,numel(x))];
        t_arcs_ct = [t_arcs_ct; 1];
        index_t = index_t + 1;
        
        % Matching condition between final computed point of backward propagated propulsed
        % arc and required final point of backward propagated propulsed arc
        Ceq = [Ceq; EquinLT_backward(1:5,end) - Equin_final_prop_back(1:5)'];
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if options.gradient == 1
            
            % Gradient constraints
            index_GCeq =  [length(x) - 13 - k : length(x) - 12 - k, length(x)-6 - k : length(x)] - (i-1)*14;
            index_GJ = [length(x) - 5 - k : length(x) - 3 - k, length(x) - k]  - (i-1)*14;
            
            eps_GCeq = eps_max;
            
            for i_GCeq = 1 : length(index_GCeq)
                
                x_var = x;
                
                x_var(index_GCeq(i_GCeq)) = x_var(index_GCeq(i_GCeq)) + delta_FD;
                
                % Azimuth, elevation
                alpha_GCeq = x_var(end-13-k-14*(i-1));
                beta_GCeq  = x_var(end-12-k-14*(i-1));
                % Initial point propagation
                start_eq_GCeq = [x_var(end-5-k-14*(i-1):end-1-k-14*(i-1)) ...
                    x_var(end-6-k-14*(i-1)) + x_var(end-k-14*(i-1))];
                % Final state
                check_eq_GCeq = x_var(end-11-k-14*(i-1):end-6-k-14*(i-1))';
                % True longitude variation
                L_prop_GCeq = linspace(start_eq_GCeq(6), check_eq_GCeq(6), n);
                
                [Ceq_var, t_var] = FABLE_grad_constr(eps_GCeq, alpha_GCeq, beta_GCeq, ...
                    L_prop_GCeq, start_eq_GCeq, check_eq_GCeq(1:5), options, constants);
                
                GCeq(temp : temp+4,index_GCeq(i_GCeq)) = ...
                    (Ceq_var - Ceq(temp:temp+4)) / delta_FD;
                
                t_arcs_var(index_t - 1, index_GCeq(i_GCeq)) = -t_var(end);
                
                % GJ for acceleration as 1/r2
                if ismember(index_GCeq(i_GCeq), index_GJ) && options.fun == 1
                    
                    B_var = sqrt(1 - x_var(end - 4 -k - (i-1)*14)^2 - ...
                        x_var(end - 3 - k- (i-1)*14)^2);
                    
                    J_var(2+2*(i-1),index_GCeq(i_GCeq)) = J_var(2+2*(i-1),index_GCeq(i_GCeq)) + ...
                        eps_max / (B_var * sqrt(constants.mu * x_var(end -5-k  - (i-1)*14) ) )* ...
                        (x_var(end  - k -  (i-1)*14));
                    
                end
                
                
            end
            
            GCeq(temp,   length(x) - 11-k - (i-1)*14) = -1;
            GCeq(temp+1, length(x) - 10-k - (i-1)*14) = -1;
            GCeq(temp+2, length(x) - 9-k - (i-1)*14) = -1;
            GCeq(temp+3, length(x) - 8-k - (i-1)*14) = -1;
            GCeq(temp+4, length(x) - 7-k - (i-1)*14) = -1;
            
            temp = temp + 5;
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if i < arcs/2
            
            % -----------------------------------------------------------------
            % Initial conditions for next forward thrust arc
            % -----------------------------------------------------------------
            Equin_initial_prop_forw = x(17+14*(i-1): 22+14*(i-1));
            Equin_final_prop_forw   = [x(23+14*(i-1): 27+14*(i-1)) ...
                x(22+14*(i-1))+x(28+14*(i-1))];
            L_prop_forw = linspace(Equin_initial_prop_forw(6), Equin_final_prop_forw(6), n);
            
            
            % -----------------------------------------------------------------
            % Initial conditions for next backward thrust arc
            % -----------------------------------------------------------------
            Equin_initial_prop_back = [x(end-19-k-14*(i-1):end-15-k-14*(i-1))...
                x(end-20-k-14*(i-1)) + x(end-14-k-14*(i-1))];
            Equin_final_prop_back   = x(end-25-k-14*(i-1):end-20-k-14*(i-1));
            L_prop_back = linspace(Equin_initial_prop_back(6), Equin_final_prop_back(6), n);
            
            
            % ---------------------------------------------------------------------
            % Backward coast arc from final point of previous low thurst arc to initial point
            % of next propulsed arc
            % ---------------------------------------------------------------------
            
            % Backward coast arc
            L_coast_back = linspace(x(end-6-k-14*(i-1)), Equin_initial_prop_back(6), n);
            Equin_initial_coast_back = x(end-11-k-14*(i-1):end-6-k-14*(i-1));
            
            % Backward propagation
            [Equin_coast_backward, t] = AnEquin_forward_m(L_coast_back, ...
                Equin_initial_coast_back, 0, 0, 0, constants.mu);
            
            % Time of flight update
            ToF = ToF - t(end);
            t_arcs = [t_arcs; -t(end)*ones(1,numel(x))];
            t_arcs_var = [t_arcs_var; -t(end)*ones(1,numel(x))];
            t_arcs_ct = [t_arcs_ct; 0];
            index_t = index_t + 1;
            
            Equin_coast = [Equin_coast, Equin_coast_backward];
            
            
            Equin_sorted{2*arcs-1-(i-1)*2} = Equin_coast_backward(:,end:-1:1);
            Equin_sorted_ct(2*arcs-1-(i-1)*2) = 0;
            
            time_sorted{2*arcs-1-(i-1)*2} = -t;
            
            if any(time_sorted{2*arcs-1-(i-1)*2})<0
                keyboard
            end
            
            % Matching condition between final point of backward coast arc and
            % initial point of next propulsed arc
            if size(Equin_coast_backward(1:5)' ) ~= size(x(end-19-k-14*(i-1):end-15-k-14*(i-1))')
                Ceq = [Ceq; Equin_coast_backward(1:5) - x(end-19-k-14*(i-1):end-15-k-14*(i-1))'];
                
            else
                Ceq = [Ceq;Equin_coast_backward(1:5)' - x(end-19-k-14*(i-1):end-15-k-14*(i-1))'];
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if options.gradient == 1
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % NEW TO COMMENT
                % Change LE
                index_GCeq = length(x) - 20-k- (i-1)*14;
                x_var = x;
                
                x_var(index_GCeq) = x_var(index_GCeq) + delta_FD;
                eps_GCeq = 0;
                alpha_GCeq = 0;
                beta_GCeq  = 0;
                
                %
                start_eq_GCeq = x_var(end-11-k-14*(i-1):end-6-k-14*(i-1));
                %
                check_eq_GCeq = [x_var(end-19-k-14*(i-1):end-15-k-14*(i-1))';...
                    x_var(end-20-k-14*(i-1))+x_var(end-14-k-14*(i-1))];
                %
                L_prop_GCeq = linspace(start_eq_GCeq(6), check_eq_GCeq(6), n);
                %
                [Ceq_var, t_var] = FABLE_grad_constr(eps_GCeq, alpha_GCeq, beta_GCeq, ...
                    L_prop_GCeq, start_eq_GCeq, check_eq_GCeq(1:5), options, constants);
                
                
                t_arcs_var(index_t - 1, index_GCeq) = -t_var(end);
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                % Gradient constraints
                index_GCeq =  [length(x) - 14-k, length(x)-11-k : length(x)-6-k] - (i-1)*14;
                eps_GCeq = 0;
                
                for i_GCeq = 1 : length(index_GCeq)
                    
                    x_var = x;
                    
                    x_var(index_GCeq(i_GCeq)) = x_var(index_GCeq(i_GCeq)) + delta_FD;
                    
                    % coast
                    alpha_GCeq = 0;
                    beta_GCeq  = 0;
                    %
                    start_eq_GCeq = x_var(end-11-k-14*(i-1):end-6-k-14*(i-1));
                    %
                    check_eq_GCeq = [x_var(end-19-k-14*(i-1):end-15-k-14*(i-1))';...
                        x_var(end-20-k-14*(i-1))+x_var(end-14-k-14*(i-1))];
                    %
                    L_prop_GCeq = linspace(start_eq_GCeq(6), check_eq_GCeq(6), n);
                    %
                    [Ceq_var, t_var] = FABLE_grad_constr(eps_GCeq, alpha_GCeq, beta_GCeq, ...
                        L_prop_GCeq, start_eq_GCeq, check_eq_GCeq(1:5), options, constants);
                    
                    GCeq(temp : temp+4,index_GCeq(i_GCeq)) = ...
                        (Ceq_var - Ceq(temp:temp+4)) / delta_FD;
                    
                    t_arcs_var(index_t - 1, index_GCeq(i_GCeq)) = -t_var(end);
                    
                end
                
                GCeq(temp,   length(x) - 19-k - (i-1)*14) = -1;
                GCeq(temp+1, length(x) - 18-k - (i-1)*14) = -1;
                GCeq(temp+2, length(x) - 17-k - (i-1)*14) = -1;
                GCeq(temp+3, length(x) - 16-k - (i-1)*14) = -1;
                GCeq(temp+4, length(x) - 15-k - (i-1)*14) = -1;
                
                
                temp = temp + 5;
                
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        end
        
        
        
    end
    
    
    if ~options.J_DeltaV
        J_out = 0;
    else
        J_out = J;
    end
    
    
    if options.gradient
        
        if options.fun == 0
            GJ = (sum(t_arcs_var(t_arcs_ct==1,:)) - J) / delta_FD;
        elseif options.fun == 1
            J_all2 = repmat(J_all, 1, numel(x));
            
            J_all2(J_var ~= 0 ) = J_var (J_var ~= 0);
            GJ = (sum(J_all2)  - J * ones(1,numel(x)) ) / delta_FD;
        end
        
    else
        GJ = [];
    end
    
    % Equality constraints
    if options.const_ToF
        
        k_ToF = 1e-3;
        Ceq = [Ceq; (ToF - curr_tof)*k_ToF];
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if options.gradient == 1 && options.const_ToF
            
            GCeq(end,:) = (sum(t_arcs_var) - sum(t_arcs)) .* k_ToF / delta_FD;
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
    end
    
    % Disequality non linear constraint
    C = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% FABLE2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif options.FABLE2
    
    % Variable to define position of variables in vector xopt;
    k = numel(x) - (8*(arcs-1)+3);
    
    % Initialise GCeq (gradient constraints)
    col_GCeq = length(x);
    rows_GCeq = 5 * arcs;
    if options.const_ToF
        GCeq = zeros(rows_GCeq + 1, col_GCeq);
    else
        GCeq = zeros(rows_GCeq, col_GCeq);
    end
    
    GJ = zeros(1, col_GCeq);
    
    delta_FD = 1e-8;
    
    t_arcs = [];
    t_arcs_var = [];
    t_arcs_ct = [];
    index_t = 1;
    
    
    count_GCeq = 1;
    
    %%
    
    % ---------------- FORWARD -------------------------------------------
    % Initial state forward propagation 
    Equin_initial_prop_forw = departure_eq;
    % Acceleration, azimuth and elevation for forward propagations
    current_eps_f = x(1); 
    current_alpha_f = x(2);
    current_beta_f = x(3);
    
    % Index of variables for GCeq of first propagated arc
    index_GCeq_f = [1 2 3];
    
    % ---------------- BACKWARD -------------------------------------------
    % Initial state backward propagation
    Equin_initial_prop_back = arrival_eq;
    % Acceleration, azimuth and elevation
    current_eps_b = x(end-9); 
    current_alpha_b = x(end-8);
    current_beta_b = x(end-7);
    
    % 
    index_GCeq_b = length(x) - [7 6 5]-k;
    
    % Angular span arcs
    DeltaL = (arrival_eq(6) - departure_eq(6) ) / (arcs);
    
    % Angular span forward and backward thrust arcs
    L_end_f = linspace(departure_eq(6), departure_eq(6) + DeltaL, n);
    L_end_b = linspace(arrival_eq(6), arrival_eq(6) - DeltaL, n);
    
    % Final state foward and backward propagation
    Equin_final_prop_forw = x(7:11);
    Equin_final_prop_back = x(end-6:end-2);
    
    
    % Pre-allocation
    J_var = zeros(arcs,numel(x));   
    J_all = [];

    for i = 1 : arcs/2
        
        %%%%%% FORWARD
        if options.fun == 0
            [EquinLT_forward, t] = AnEquin_forward_m(L_end_f, ...
                Equin_initial_prop_forw', current_eps_f, current_alpha_f, ...
                current_beta_f, constants.mu);
        elseif options.fun == 1
            [EquinLT_forward, t] = AnEquin_forward_m_r2(L_end_f, ...
                Equin_initial_prop_forw', current_eps_f, current_alpha_f, ...
                current_beta_f, constants.mu);
        end
        
        
        Equin_thrust = [Equin_thrust, EquinLT_forward];
        if any(Equin_thrust(1,:)<10)
            keyboard
        end
        
        Equin_sorted{i} = EquinLT_forward;
        Equin_sorted_ct(i) = 1;
        
        time_sorted{i} = t;
        
        if options.fun == 0
            % Cost function: time spent using propulsion
            J = J + t(end) * current_eps_f;
        elseif options.fun == 1
            % Cost function: analytic expression DeltaV
            B = sqrt(1 - Equin_initial_prop_forw(2)^2 - Equin_initial_prop_forw(3)^2);
            J = J + current_eps_f / (B * sqrt(constants.mu * Equin_initial_prop_forw(1)) )* ...
                (L_end_f(end) - L_end_f(1));
            J_all = [J_all; current_eps_f / (B * sqrt(constants.mu * Equin_initial_prop_forw(1)) )* ...
                (L_end_f(end) - L_end_f(1))];
        end
        
        % Time of flight update
        ToF = ToF + t(end);
        
        t_arcs = [t_arcs; t(end)*ones(1,numel(x))];
        t_arcs_var = [t_arcs_var; t(end)*ones(1,numel(x))];
        t_arcs_ct = [t_arcs_ct; 0];
        
        index_t = index_t + 1;
        
%         T_arcs = t(end);
        
        
        % Matching condition between final point of forward propagated propulsed
        % arc and required final point of forward propagated propulsed arc
        Ceq = [Ceq; EquinLT_forward(1:5,end) - Equin_final_prop_forw(1:5)'];
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if options.gradient == 1
            
            if i == 1
                index_GJ = 1;
            else
                index_GJ = [4, 7, 8, 9] + (i-2) * 8;
            end
            
            if i == 1 && numel(x) == 8*(arcs-1)+3 + 2
                % Change alpha and get effect on GCeq
                % Azimuth and declination
                alpha_var = x(end-1) + delta_FD;
                delta_var = x(end) + delta_FD;
                
                % ------------alpha variation
                departure_vel_var = parameters.dep_body_v + ...
                    v_inf * [cos(alpha_var) * cos(delta), sin(alpha_var) * cos(delta), sin(delta)];
                
                departure_kep_var = cart2kep([departure_pos, departure_vel_var], constants.mu);
                
                departure_eq_var = kep2eq(departure_kep_var);

                % Thrust arc
                eps_GCeq = current_eps_f;
                alpha_GCeq = current_alpha_f;
                beta_GCeq = current_beta_f;
                % Final point propagation
                L_prop_GCeq = L_end_f;
                % Initial point propagation
                start_eq_GCeq = departure_eq_var;
                % Final state of propagation should be...:
                check_eq_GCeq = Equin_final_prop_forw(1:5);
                % Constraint under x_var
                [Ceq_var, t_var] = FABLE_grad_constr(eps_GCeq, alpha_GCeq, beta_GCeq, L_prop_GCeq, ...
                    start_eq_GCeq, check_eq_GCeq', options, constants);
                % GCeq
                GCeq(1:5,end-1) = (Ceq_var - Ceq(1:5)) / delta_FD;
                
                %
                t_arcs_var(index_t - 1, end-1) = t_var(end);
                
                
                % Variation due to alpha on the objective
                B_var = sqrt(1 - departure_eq_var(2)^2 - ...
                    departure_eq_var(3)^2);
                
                J_var(1,end-1) = ...
                    J_var(1,end-1) + ...
                    eps_GCeq /....
                    (B_var * sqrt(constants.mu * departure_eq_var(1)) )* ...
                    (L_end_f(end) - L_end_f(1));
                
                
                % -----------delta variation
                departure_vel_var = parameters.dep_body_v + ...
                    v_inf * [cos(alpha) * cos(delta_var), sin(alpha) * cos(delta_var), sin(delta_var)];
                
                departure_kep_var = cart2kep([departure_pos, departure_vel_var], constants.mu);
                
                departure_eq_var = kep2eq(departure_kep_var);
                
                % Thrust arc
                eps_GCeq = current_eps_f;
                alpha_GCeq = current_alpha_f;
                beta_GCeq = current_beta_f;
                % Final point propagation
                L_prop_GCeq = L_end_f;
                % Initial point propagation
                start_eq_GCeq = departure_eq_var;
                % Final state of propagation should be...:
                check_eq_GCeq = Equin_final_prop_forw(1:5);
                % Constraint under x_var
                [Ceq_var, t_var] = FABLE_grad_constr(eps_GCeq, alpha_GCeq, beta_GCeq, L_prop_GCeq, ...
                    start_eq_GCeq, check_eq_GCeq', options, constants);
                % GCeq
                GCeq(1:5,end) = (Ceq_var - Ceq(1:5)) / delta_FD;
                
                t_arcs_var(index_t - 1, end) = t_var(end);
                
                % Variation due to delta on the objective
                B_var = sqrt(1 - departure_eq_var(2)^2 - ...
                    departure_eq_var(3)^2);
                
                J_var(1,end) = ...
                    J_var(1,end) + ...
                    eps_GCeq /....
                    (B_var * sqrt(constants.mu * departure_eq_var(1)) )* ...
                    (L_end_f(end) - L_end_f(1));
            end
            
            for i_GCeq = 1 : length(index_GCeq_f)
                
                x_var = x;
                x_var(index_GCeq_f(i_GCeq)) = x_var(index_GCeq_f(i_GCeq)) + ...
                    delta_FD;
                
                if i == 1
                    eps_GCeq = x_var(1); 
                    alpha_GCeq = x_var(2);
                    beta_GCeq  = x_var(3);
                    start_eq_GCeq = departure_eq;
%                     if index_GCeq_f(i_GCeq) == 1
%                         GJ(index_GCeq_f(i_GCeq)) = GJ(index_GCeq_f(i_GCeq)) + t(end);
%                     end
                else
                    eps_GCeq = x_var(4 + (i-2)*8); 
                    alpha_GCeq = x_var(5 + (i-2)*8);
                    beta_GCeq  = x_var(6 + (i-2)*8);
                    start_eq_GCeq = [x_var(7+8*(i-2): 11+8*(i-2)) L_end_f(1)];
%                     if index_GCeq_f(i_GCeq) == 4 + (i-2)*8
%                         GJ(index_GCeq_f(i_GCeq)) = GJ(index_GCeq_f(i_GCeq)) + t(end);
%                     end
                end

                check_eq_GCeq = x_var(7+8*(i-1): 11+8*(i-1));
                
                L_prop_GCeq = L_end_f;
                
                [Ceq_var, t_var] = FABLE_grad_constr(eps_GCeq, alpha_GCeq, beta_GCeq,...
                    L_prop_GCeq, start_eq_GCeq, check_eq_GCeq(1:5)', options, constants);
                
                GCeq(count_GCeq:count_GCeq+4,index_GCeq_f(i_GCeq)) = ...
                    (Ceq_var - Ceq(count_GCeq:count_GCeq+4)) / delta_FD;
                
                t_arcs_var(index_t - 1, index_GCeq_f(i_GCeq)) = t_var(end);
            
                % GJ for acceleration as 1/r2
                if ismember(index_GCeq_f(i_GCeq), index_GJ) && options.fun == 1
                    
                    if i == 1
                        
                        % Variation due to epsilon 1
                        B = sqrt(1 - departure_eq(2)^2 - ...
                            departure_eq(3)^2);
                        
                        J_var(1,1) = ...
                            J_var(1,1) + ...
                            x_var(1 ) /....
                            (B * sqrt(constants.mu * departure_eq(1)) )* ...
                            (L_end_f(end) - L_end_f(1));
                        
                    else

                        B_var = sqrt(1 - x_var(8 + (i-2)*8)^2 - ...
                            x_var(9 + (i-2)*8)^2);
                        
                        J_var(1+2*(i-1),index_GCeq_f(i_GCeq)) = ...
                            J_var(1+2*(i-1),index_GCeq_f(i_GCeq)) + ...
                            x_var(4 + (i-2)*8) /....
                            (B_var * sqrt(constants.mu * x_var(7 + (i-2)*8)) )* ...
                            (L_end_f(end) - L_end_f(1));
                    end
                end
                %                 GJ(index_GCeq_f(i_GCeq)) = GJ(index_GCeq_f(i_GCeq)) + ...
%                     eps_GCeq * ( time_var(end) - t(end) )/ delta_FD;
                
            end
            
            GCeq(count_GCeq,     7 + (i-1)*8) = -1;
            GCeq(count_GCeq + 1, 8 + (i-1)*8) = -1;
            GCeq(count_GCeq + 2, 9 + (i-1)*8) = -1;
            GCeq(count_GCeq + 3, 10 + (i-1)*8) = -1;
            GCeq(count_GCeq + 4, 11 + (i-1)*8) = -1;
            
            count_GCeq = count_GCeq + 5;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Equin_initial_prop_forw = [x(7 + (i-1)*8 : 11 + (i-1)*8 ) L_end_f(end)];
        Equin_final_prop_forw = x(15 + (i-1)*8 : 19 + (i-1)*8 );
        current_eps_f = x(4 + (i-1)*8); %* (parameters.UB(1) - parameters.LB(1)) + parameters.LB(1);
        current_alpha_f =  x(5 + (i-1)*8);
        current_beta_f = x(6 + (i-1)*8);
        L_end_f = linspace(L_end_f(end), L_end_f(end) + DeltaL, n);
        index_GCeq_f = [4:1:11] + (i-1)*8;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        % ---------------------------------------------------------------------
        % Backward propagated low-thrust arc
        % ---------------------------------------------------------------------
        
        if options.fun == 0
            [EquinLT_backward, t] =  AnEquin_forward_m(L_end_b, ...
                Equin_initial_prop_back', current_eps_b,...
                current_alpha_b, current_beta_b, ...
                constants.mu);
        elseif options.fun == 1
            [EquinLT_backward, t] =  AnEquin_forward_m_r2(L_end_b, ...
                Equin_initial_prop_back', current_eps_b,...
                current_alpha_b, current_beta_b, ...
                constants.mu);
        end

        Equin_thrust = [Equin_thrust, EquinLT_backward];
        if any(Equin_thrust(1,:)<10)
            keyboard
        end
        
        Equin_sorted{arcs-i+1} = EquinLT_backward(:,end:-1:1);
        Equin_sorted_ct(arcs-i+1) = 1;
        
        time_sorted{arcs-i+1} = -t;
        
        if options.fun == 0
            % Cost function: time spent using propulsion
                J = J - t(end) * current_eps_b;
        elseif options.fun == 1
            % Cost function: analytic computation DeltaV
            B = sqrt(1 - Equin_initial_prop_back(2)^2 - Equin_initial_prop_back(3)^2);
            J = J + current_eps_b / (B * sqrt(constants.mu * Equin_initial_prop_back(1)) )* ...
                (L_end_b(1) - L_end_b(end));
            J_all = [J_all; current_eps_b / (B * sqrt(constants.mu * Equin_initial_prop_back(1)) )* ...
                (L_end_b(1) - L_end_b(end))];
        end
        
        
        
        % Time of flight update
        ToF = ToF - t(end);
        
        t_arcs = [t_arcs; -t(end)*ones(1,numel(x))];
        t_arcs_var = [t_arcs_var; -t(end)*ones(1,numel(x))];
        t_arcs_ct = [t_arcs_ct; 1];
        index_t = index_t + 1;
        
        % Matching condition between final computed point of backward propagated propulsed
        % arc and required final point of backward propagated propulsed arc
        Ceq = [Ceq; EquinLT_backward(1:5,end) - Equin_final_prop_back(1:5)'];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if options.gradient == 1
            
            if i == 1
                index_GJ = numel(x) - k - 7;
            else
                a_temp = numel(x) - k - 15 - (i-2)*8;
                b_temp = numel(x) - k -4 -(i-2)*8: numel(x) - k-(i-2)*8;
                index_GJ = [a_temp, b_temp];
            end
            
            
            % Gradient constraints
            
            for i_GCeq = 1 : length(index_GCeq_b)
                x_var = x;
                x_var(index_GCeq_b(i_GCeq)) = x_var(index_GCeq_b(i_GCeq)) + ...
                    delta_FD;
                
                
                eps_GCeq   = x_var(end - k - 7 - (i-1)*8); 
                alpha_GCeq = x_var(end - k - 6 - (i-1)*8);
                beta_GCeq  = x_var(end - k - 5 - (i-1)*8);
                
                if i == 1
                    start_eq_GCeq = arrival_eq;
                else
                    start_eq_GCeq = [x_var(end - 4 - k - 8*(i-2) : ...
                        end - k -8*(i-2)) L_end_b(1)];
                end

      
                check_eq_GCeq = x_var(end - 4 - k - 8*(i-1): end - k - 8*(i-1));
                
                L_prop_GCeq = L_end_b;
                
                [Ceq_var, t_var] = FABLE_grad_constr(eps_GCeq, alpha_GCeq, beta_GCeq,...
                    L_prop_GCeq, start_eq_GCeq, check_eq_GCeq(1:5)', options, constants);
                
                GCeq(count_GCeq:count_GCeq+4,index_GCeq_b(i_GCeq)) = ...
                    (Ceq_var - Ceq(count_GCeq:count_GCeq+4)) / delta_FD;
                
                %                 GJ(index_GCeq_b(i_GCeq)) = GJ(index_GCeq_b(i_GCeq)) + ...
                %                     eps_GCeq * ( -time_var(end) + t(end) )/ delta;
                
                t_arcs_var(index_t - 1, index_GCeq_b(i_GCeq)) = -t_var(end);
                
                % GJ for acceleration as 1/r2
                if ismember(index_GCeq_b(i_GCeq), index_GJ) && options.fun == 1
                    
                    if i == 1
                        B_var = sqrt(1-arrival_eq(2)^2 - arrival_eq(3)^2);
                        J_var(2+2*(i-1),index_GCeq_b(i_GCeq)) =...
                            J_var(2+2*(i-1),index_GCeq_b(i_GCeq)) + ...
                            x_var(end - k - 7) / ...
                            (B_var * sqrt(constants.mu * arrival_eq(1) ) )* ...
                            (L_end_b(1) - L_end_b(end));
                    else
                        
                        B_var = sqrt(1 - x_var(end - 3 -k - (i-2)*8)^2 - ...
                            x_var(end - 2 - k- (i-2)*8)^2);
                        
                        J_var(2+2*(i-1),index_GCeq_b(i_GCeq)) =...
                            J_var(2+2*(i-1),index_GCeq_b(i_GCeq)) + ...
                            x_var(end -7 - k - (i-1)*8)/ ...
                            (B_var * sqrt(constants.mu * x_var(end -4-k  - (i-2)*8) ) )* ...
                            (L_end_b(1) - L_end_b(end));
                    end
                end
                
                
            end
            GCeq(count_GCeq,   end - k - 4 - (i-1)*8) = -1;
            GCeq(count_GCeq+1, end - k - 3 - (i-1)*8) = -1;
            GCeq(count_GCeq+2, end - k - 2 - (i-1)*8) = -1;
            GCeq(count_GCeq+3, end - k - 1 - (i-1)*8) = -1;
            GCeq(count_GCeq+4, end - k -     (i-1)*8) = -1;
            
            %         GJ(index_GCeq_f(i_GCeq)) =  GJ(index_GCeq_f(i_GCeq)) + ...
            %             current_eps_b * ( time_var(end) - t(end) )/ delta ;
            
            count_GCeq = count_GCeq + 5;
            
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Initial condition for propagation of next backward arc
        Equin_initial_prop_back = [x(end - 4 - k - (i-1)*8 : ...
            end - k - (i-1)*8 ) L_end_b(end)];
        % Final desired state for backward propagation of next arc
        Equin_final_prop_back = x(end - 12 - k - (i-1)*8 : ...
            end - 8 - k - (i-1)*8 );
        % Acceleration, azimuth and elevation for next arc
        current_eps_b   = x(end - 15 - k - (i-1)*8); 
        current_alpha_b = x(end - 14 - k - (i-1)*8);
        current_beta_b  = x(end - 13 - k - (i-1)*8);
        % Range of L for next backward propagated arc
        L_end_b = linspace(L_end_b(end), L_end_b(end) - DeltaL, n);
        % Index for GCeq for next arc
        a_temp = length(x) - 15 - k - (i-1) * 8 : length(x) - 13 - k - (i-1) * 8;
        b_temp = length(x) - 4 - k - (i-1) * 8 : length(x)  - k -(i-1) * 8;
        index_GCeq_b = [a_temp, b_temp];
        
        
    end
    
    
    
    if ~options.J_DeltaV
        J_out = 0;
    else
        J_out = J;
        %     J_out = ToF;
    end
    
    
    if options.gradient
        
        if options.fun == 0
            GJ = (sum(t_arcs_var(t_arcs_ct==1,:)) - J) / delta_FD;
        elseif options.fun == 1
            J_all2 = repmat(J_all, 1, numel(x));
            
            J_all2(J_var ~= 0 ) = J_var (J_var ~= 0);
            GJ = (sum(J_all2)  - J * ones(1,numel(x)) ) / delta_FD;
        end
        
    else
        GJ = [];
    end
    
    
    % Equality constraints
    if options.const_ToF
        
        k_ToF = 1e-3;
        Ceq = [Ceq; (ToF - curr_tof)*k_ToF];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if options.gradient == 1 && options.const_ToF
            
            GCeq(end,:) = (sum(t_arcs_var) - sum(t_arcs)) .* k_ToF / delta_FD;
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
    end
    
    % Disequality non linear constraint
    C = [];
    
    
end

GC = [];


%% Plot

if parameters.plot_flag
    
    % Thrusting and coasting arcs
    Cart_thrust = zeros(length(Equin_thrust),6);
    Cart_coast  = zeros(length(Equin_coast),6);
    for k = 1 : length(Equin_coast)
        Cart_coast(k,:) = eq2cart_m(Equin_coast(:,k),constants.mu);
    end
    for k = 1 : length(Equin_thrust)
        Cart_thrust(k,:) = eq2cart_m(Equin_thrust(:,k),constants.mu);
    end
    
    % Departure orbit
    departure_cart = eq2cart_m(departure_eq, constants.mu);
    
    % Arrival orbit
    arrival_cart = eq2cart_m(arrival_eq, constants.mu);
    
    
    % Period of departure and arrival trajectories
    T_departure = 2 * pi * sqrt(departure_eq(1)^3/constants.mu);
    T_arrival = 2 * pi * sqrt(arrival_eq(1)^3/constants.mu);
    
    figure
    hold on
    % Plot departure and arrival positions
    plot3(departure_cart(1)*constants.DU/constants.AU,...
        departure_cart(2)*constants.DU/constants.AU,...
        departure_cart(3)*constants.DU/constants.AU,...
        'bo','MarkerFaceColor','b','MarkerSize',8)
    plot3(arrival_cart(1)*constants.DU/constants.AU,...
        arrival_cart(2)*constants.DU/constants.AU,...
        arrival_cart(3)*constants.DU/constants.AU,...
        'ko','MarkerFaceColor','k','MarkerSize',8)
    
    % Plot departure and arrival trajectories
    plot_trajectory(departure_cart(1:3)*constants.DU/constants.AU,...
        departure_cart(4:6)*constants.DU/constants.AU/(constants.TU/constants.TU_AU),...
        T_departure*constants.TU/constants.TU_AU,constants.mu,'b',2);
    plot_trajectory(arrival_cart(1:3)*constants.DU/constants.AU,...
        arrival_cart(4:6)*constants.DU/constants.AU/(constants.TU/constants.TU_AU),...
        T_arrival*constants.TU/constants.TU_AU,constants.mu,'k',2);
    
    % Plot thrust and coast arcs
    index = 1;
    while index <= length(Equin_coast)
        plot3(Cart_coast(index : index + parameters.n - 2, 1)*constants.DU/constants.AU,...
            Cart_coast(index : index + parameters.n - 2,2)*constants.DU/constants.AU, ...
            Cart_coast(index : index + parameters.n - 2,3)*constants.DU/constants.AU,'g-','LineWidth',2);
        index = index + parameters.n-1;
    end
    
    index = 1;
    while index <= length(Equin_thrust)
        plot3(Cart_thrust(index : index + parameters.n - 2, 1)*constants.DU/constants.AU,...
            Cart_thrust(index : index + parameters.n - 2,2)*constants.DU/constants.AU, ...
            Cart_thrust(index : index + parameters.n - 2,3)*constants.DU/constants.AU,'r-','LineWidth',2);
        index = index + parameters.n-1;
    end
    grid on
    xlabel('x [AU]')
    ylabel('y [AU]')
    zlabel('z [AU]')
    view(0,90)
    axis equal
    legend('Departure Point','Arrival Point','Initial orbit','Final orbit')
    
    if options.FABLE2
        
        ls = 14;
        
        figure
        hold on
        index = 1;
        count = 1;
        eps_value = x(1) / (constants.TU*constants.sec_day)^2 * constants.DU*1000;
        while count <= length(Equin_sorted)
            plot(Equin_sorted{count}(6,:)*180/pi,...
                eps_value * ones(1, length(Equin_sorted{count}(6,:))), ...
                'LineWidth',2);
            eps_value = x(4 + (count-1) * 8)/ (constants.TU*constants.sec_day)^2 * constants.DU*1000;
            count = count + 1;
        end
        grid on
        xlabel('L [deg]','FontSize',ls)
        ylabel('\epsilon [m/s^2]','FontSize',ls)
        set(gca,'FontSize',ls)
        
        
        
        figure
        hold on
        count = 1;
        beta_value = x(3) *180/pi;
        plot(Equin_sorted{count}(6,:)*180/pi,...
            beta_value * ones(1, length(Equin_sorted{count*2}(6,:))), ...
            'LineWidth',2);
        while count < length(Equin_sorted)
            
            beta_value = x(6 + (count-1) * 8) * 180/pi;
            count = count + 1;
            plot(Equin_sorted{count}(6,:)*180/pi,...
                beta_value * ones(1, length(Equin_sorted{count}(6,:))), ...
                'LineWidth',2);
        end
        grid on
        xlabel('L [deg]','FontSize',ls)
        ylabel('\beta [deg]','FontSize',ls)
        set(gca,'FontSize',ls)
        

        
        figure
        hold on
        count = 1;
        alpha_value = x(2) *180/pi;
        while count <= length(Equin_sorted)
            plot(Equin_sorted{count}(6,:)*180/pi,...
                alpha_value * ones(1, length(Equin_sorted{count}(6,:))), ...
                'LineWidth',2);
            alpha_value = x(5 + (count-1) * 8) * 180/pi;
            count = count + 1;
        end
        grid on
        xlabel('L [deg]','FontSize',ls)
        ylabel('\alpha [deg]','FontSize',ls)
        set(gca,'FontSize',ls)
        
        
    end
    
end






end
