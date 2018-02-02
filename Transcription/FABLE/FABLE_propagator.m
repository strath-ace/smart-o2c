function [ToF, x_final, J, output] = ...
    FABLE_propagator(L_f, x, u, parameters, options, constants)

% =========================================================================
% Forward propagation using analytic equations of a solution vector
% obtained using FABLE. 
% =========================================================================

% Input: L_f        -> final true longitude
%        x          -> 7x1 vector of initial conditions(6 equinoctial elements, mass)
%        u          -> matrix of controls. Number of rows is equal to number of thrust
%                   arcs. For FABLE1: 1st column: true longitude ON node. 2nd column: true
%                   longitude OFF node. 3rd column: azimuth angle. 4th columns:
%                   elevation angle.
%                   For FABLE2: 1st column: acceleration. 2nd column:
%                   azimuth angle. 3rd column: elevation angle.
%        parameters -> structure containing the parameters of the problem
%        options    -> structure containing the options for the problem
%        constants  -> structure of constants
% 
% Output: ToF     -> time of flight of the transfer
%         x_final -> final equinoctial elements
%         J       -> objective function for FABLE SS
%         output  -> structure containing orbital elements variation
% 
% Marilena Di Carlo, marilena.di-carlo@strath.ac.uk

%% Initialisation

% Mass
m = x(7);

% Number of transfer arcs
arcs         = parameters.arcs;

% Number of steps for analytic integration
n            = parameters.n;

% Departure state
departure_eq = x(1:6);

% Cost function
J = 0;

% Time of flight
ToF = 0;

% Initialize varaibles
Equin_coast = [];
Equin_thrust = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FABLE 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if options.FABLE1
    
    % Low-thrust acceleration
    eps      = parameters.T_adim / parameters.m;
    
    %% FIRST PROPULSED ARC
    
    % Initial equinoctial element for forward propagation on the first low
    % thrust arc
    L_initial_prop_forw  = u(1,1);
    
    % Final equinoctial elements of first forward propagation on the first low
    % thrust arc
    L_final_prop_forw = u(1,2);
    

    %% FIRST COAST ARC?
    
    % Define an initial coast arc only if the first point of the vector x
    % does not coincide with the departure point
    
    if departure_eq(6) == L_initial_prop_forw
        
        % If the longitude of the initial point of the first LT arc coincides with the
        % departure longitude, there is no need of using a coast arc and the
        % first equality contraint (coincidence of the two points) has to be
        % defined
        
        % True longitude variation over first low thrust arc
        L_prop_forw  = linspace(L_initial_prop_forw,  L_final_prop_forw, n);
        
        Equin_initial_prop_forw  = departure_eq;
        
    else
        % If the initial point of the first propulsed arc does not coincide
        % with the departure point, use a coast arc to reach the first point
        
        
        % True longitude variation
        L_coast_forw = linspace(departure_eq(6), L_initial_prop_forw, n);
        
        
        % Propagation. Since this is a coast arc it is not necessary to choose
        % the model for the thrust. Therefore the following is commented and
        % substituted by
        [Equin_coast_forw, t] = AnEquin_forward_m(L_coast_forw, departure_eq, ...
            0, 0, 0, constants.mu);
        
        Equin_coast = [Equin_coast, Equin_coast_forw];
        
        Equin_sorted{1} = Equin_coast_forw;
        Equin_sorted{1} = [Equin_sorted{1}; ...
            m * ones(1,size(Equin_sorted{1},2))];
        
        Equin_sorted_ct(1) = 0;
        
        time_sorted{1} = t;
        
        
        % Time of flight update
        ToF = ToF + t(end);
        
        % The final point of the initial coast arc has to coincide with the initial point
        % of the first low thrust arc
        
        Equin_initial_prop_forw = Equin_coast_forw(:,end);
        
        L_prop_forw  = linspace(Equin_initial_prop_forw(6),  L_final_prop_forw, n);
        
    end
    
    %% For cycle over number of thrust arcs
    
    for i = 1 : arcs

        % ---------------------------------------------------------------------
        % Forward propagated low-thrust arc
        % ---------------------------------------------------------------------
        if options.fun == 0
            [EquinLT_forward, t] = AnEquin_forward_m(L_prop_forw,...
                Equin_initial_prop_forw', eps, u(i,3), u(i,4), ...
                constants.mu);
        elseif options.fun == 1
            [EquinLT_forward, t] = AnEquin_forward_m_r2(L_prop_forw,...
                Equin_initial_prop_forw', eps, u(i,3), u(i,4), ...
                constants.mu);
        end
        
        
        if i < arcs
            % Next Forward coast arc
            Equin_initial_coast_forw = EquinLT_forward(:,end);
            L_coast_forw = linspace(Equin_initial_coast_forw(6), u(i+1,1), n);
        else
            % Next Forward coast arc
            Equin_initial_coast_forw = EquinLT_forward(:,end);
            L_coast_forw = linspace(Equin_initial_coast_forw(6), L_f, n);
        end
        
        Equin_thrust = [Equin_thrust, EquinLT_forward];
        
        Equin_sorted{2+(i-1)*2} = EquinLT_forward;
        Equin_sorted_ct(2+(i-1)*2) = 1;
        
        time_sorted{2+(i-1)*2} = t;
        
        if options.fun == 0
            % Cost function: time spent using propulsion
            J = J + t(end);
            % Update mass
            kep_temp = eq2kep(Equin_initial_prop_forw);
            cart_temp = kep2cart(kep_temp, constants.mu);
            r_temp = norm(cart_temp(1:3));
            m_temp = m - parameters.T_adim / r_temp^2 / (parameters.Isp_adim * constants.g0_adim) * t;
            m = m_temp(end);
            Equin_sorted{2+(i-1)*2} = [Equin_sorted{2+(i-1)*2}; m_temp];
        elseif options.fun == 1
            % Cost function: analytic expression DeltaV
            B = sqrt(1 - Equin_initial_prop_forw(2)^2 - Equin_initial_prop_forw(3)^2);
%             J = J + eps / (B * sqrt(constants.mu * Equin_initial_prop_forw(1)) )* ...
%                 (L_prop_forw(end) - L_prop_forw(1));
            J = J + parameters.T_adim/m / (B * sqrt(constants.mu * Equin_initial_prop_forw(1)) )* ...
                (L_prop_forw(end) - L_prop_forw(1));
            % Update mass
            %         m =m- m * eps / (B * sqrt(constants.mu * Equin_initial_prop_forw(1)) )* ...
            %             (L_prop_forw(end) - L_prop_forw(1)) / (parameters.Isp_adim * constants.g0_adim);
            m_temp =m- parameters.T_adim / (B * sqrt(constants.mu * Equin_initial_prop_forw(1)) )* ...
                (L_prop_forw - L_prop_forw(1)) / (parameters.Isp_adim * constants.g0_adim);
            m = m_temp(end);
            Equin_sorted{2+(i-1)*2} = [Equin_sorted{2+(i-1)*2}; m_temp(2:end)];
        end
        
        
        
        % Update acceleration
        eps = parameters.T_adim / m;
        
        % Time of flight update
        ToF = ToF + t(end);
        
        % ---------------------------------------------------------------------
        % Forward coast arc from final point of low thurst arc to initial point
        % of next propulsed arc
        % ---------------------------------------------------------------------
        [Equin_coast_forward, t] = AnEquin_forward_m(L_coast_forw, ...
            Equin_initial_coast_forw, 0, 0, 0, constants.mu);
        
        
        
        if i < arcs
            % Next thrust arc
            Equin_initial_prop_forw = Equin_coast_forward(:,end);
            L_prop_forw  = linspace(Equin_initial_prop_forw(6),  ...
                u(i+1,2), n);
            
        end
        
        Equin_coast = [Equin_coast, Equin_coast_forward];
        
        
        Equin_sorted{3+(i-1)*2} = Equin_coast_forward;
        Equin_sorted{3+(i-1)*2} = [Equin_sorted{3+(i-1)*2}; ...
            m * ones(1,size(Equin_sorted{3+(i-1)*2},2))];
        
        Equin_sorted_ct(3+(i-1)*2) = 0;
        
        time_sorted{3+(i-1)*2} = t;
        
        if any(time_sorted{3+(i-1)*2})<0
            keyboard
        end
        
        % Time of flight update
        ToF = ToF + t(end);
        
        
        
        
    end
    
    
    x_end = Equin_coast(:,end);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FABLE 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif options.FABLE2
    
%     arrival_eq   = parameters.arrival_eq;
    
    DeltaL = (L_f - departure_eq(6) ) / (arcs);
   
    Equin_initial_prop = departure_eq;
    L_end = linspace(departure_eq(6), departure_eq(6) + DeltaL, n);
    

    %%
    
    Equin_initial_prop_forw = departure_eq;
    
    for i = 1 : arcs
        
        % Control
        current_eps   = u(i,1);
        current_alpha = u(i,2);
        current_beta  = u(i,3);
        
 
        if options.fun == 0
            [EquinLT_forward, t] = AnEquin_forward_m(L_end, ...
                Equin_initial_prop', current_eps, current_alpha, ...
                current_beta, constants.mu);
        elseif options.fun == 1
            [EquinLT_forward, t] = AnEquin_forward_m_r2(L_end, ...
                Equin_initial_prop', current_eps, current_alpha, ...
                current_beta, constants.mu);
        end
        
        
        Equin_thrust = [Equin_thrust, EquinLT_forward];
        
        Equin_sorted{i} = EquinLT_forward;
        Equin_sorted_ct(i) = 1;
        
        time_sorted{i} = t;
        
        % Cost function: time spent using propulsion
        if options.fun == 0
            J = J + t(end) * current_eps;
            % Update mass
            kep_temp = eq2kep(Equin_initial_prop_forw);
            cart_temp = kep2cart(kep_temp, constants.mu);
            r_temp = norm(cart_temp(1:3));
            m = m - parameters.T_adim / r_temp^2 / (parameters.Isp_adim * constants.g0_adim) * t(end);
        elseif options.fun == 1
            % Cost function: analytic expression DeltaV
            B = sqrt(1 - Equin_initial_prop(2)^2 - Equin_initial_prop(3)^2);
            J = J + current_eps / (B * sqrt(constants.mu * Equin_initial_prop(1)) )* ...
                (L_end(end) - L_end(1));
            % Update mass
            m_temp = m- m * current_eps / (B * sqrt(constants.mu * Equin_initial_prop(1)) )* ...
                (L_end - L_end(1)) / (parameters.Isp_adim * constants.g0_adim);
            m = m_temp(end);
        end
        
        % Time of flight update
        ToF = ToF + t(end);
        
        
        % Initial conditions next arc
        Equin_initial_prop = EquinLT_forward(:,end);
        L_end = linspace(L_end(end), L_end(end) + DeltaL, n);
        
        
        
    end
    
    x_end = EquinLT_forward(:,end);
end

%% Outputs common to FABLE 1 and FABLE2

if ~options.J_DeltaV
    J_out = 0;
else
    J_out = J;
end


x_final = [x_end; m];


output.Equin_sorted = Equin_sorted;
output.Equin_sorted_ct = Equin_sorted_ct;
output.time_sorted = time_sorted;

output.Equin_thrust = Equin_thrust;
output.Equin_coast = Equin_coast;
