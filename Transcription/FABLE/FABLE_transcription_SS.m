function [J_out,C,Ceq, GJ, GC, GCeq, ToF,Equin_sorted,Equin_sorted_ct,time_sorted,J,E0] = ...
    FABLE_transcription_SS(x, parameters, options, constants, inputs)

% =========================================================================
% Single shooting transcription of the low-thrust problem for FABLE.
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

% The vector x is 1x4*arcs. 
% For each arc, 4 variable are defined: 
% - alpha (azimuth)
% - beta(elevation)
% - true longitude of the initial point of the propulsed arc 
% - true longitude of the final point of the propulsed arc.

% They are organized in the x vector as follows:
% [alpha, beta, [L], [Delta L]]
% for each propulsed arc. 
% =========================================================================
% Author: Marilena Di Carlo
% email: marilena.di-carlo@strath.ac.uk
% =========================================================================



%% Initialisation

% Number of transfer arcs
arcs         = parameters.arcs;


% Equinoctial elements arrival point: 
% 1st option: the time of flight is fixed, arrival point is defined
if options.const_ToF == 1
    arrival_eq   = parameters.arrival_eq;
else
    % Compute arrival position based on time of fligth
    Apophis = inputs.Apophis;
    % 2nd option: recompute arrival position based on ToF
    % Adimensional semimajor axis    
    kk = fix( ( (inputs.departure_date_MJD2000 + x(end)*constants.TU - Apophis.epochMJD2000 ) / constants.TU) ...
        /  Apophis.T );
    
    % Mean anomaly at arrival
    Apophis.M = mod(Apophis.M0 + ...
        Apophis.n * (x(end)*constants.TU + inputs.departure_date_MJD2000- Apophis.epochMJD2000)/constants.TU + ...
        -2 * kk * pi, 2*pi);
    
    % Eccentric anomaly at arrival
    Apophis.E = kepler_E(Apophis.e, Apophis.M);
    
    % True anomaly at arrival
    Apophis.cos_theta = ( cos(Apophis.E) - Apophis.e ) / ...
        (1 - Apophis.e * cos(Apophis.E));
    
    Apophis.sin_theta = ( sin(Apophis.E)  * sqrt(1 - Apophis.e^2) ) / ...
        (1 - Apophis.e * cos(Apophis.E));
    
    Apophis.theta = atan2(Apophis.sin_theta,Apophis.cos_theta);
    
    % Orbital element of apophis at arrival
    arrival_kep = [Apophis.a, Apophis.e, Apophis.i, Apophis.RAAN, Apophis.omega, ...
        Apophis.theta];
    
    arrival_eq = kep2eq(arrival_kep);
    
    % Increase the arrival true longitude by 2*pi times the number of
    % revolutions
    arrival_eq(6) = arrival_eq(6) + 2 * pi * (options.n_rev);
end




% Variable to define position of variables in vector xopt;
k1 = options.const_ToF ==  0;



% Time of flight
curr_tof     = parameters.curr_tof;


% Equinoctial elements departure point
if (options.FABLE1 && numel(x) == 4*arcs && options.const_ToF == 1) || ...
   (options.FABLE1 && numel(x) == 4*arcs+1 && options.const_ToF == 0) || ...
    (options.FABLE2 && numel(x) == 3*arcs)
    % If departure state is given and declination and azimuth at departure
    % do not have to be optimised
    departure_eq = parameters.departure_eq;
elseif (options.FABLE1 && numel(x) == 4*arcs + 2 && options.const_ToF == 1) || ...
       (options.FABLE1 && numel(x) == 4*arcs + 3 && options.const_ToF == 0) || ...
        options.FABLE2 && numel(x) == 3*arcs + 2
    % Departure position
    departure_pos = parameters.dep_body_r;
    % Velocity at departure is given by the sum of the velocity of the
    % central body and the velocity vinf
    v_inf = parameters.vinf0/ constants.DU * constants.TU * constants.sec_day;
    % Azimuth and declination
    alpha = x(end-k1-1);
    delta = x(end-k1);
    
    departure_vel = parameters.dep_body_v + ...
        v_inf * [cos(alpha) * cos(delta), sin(alpha) * cos(delta), sin(delta)];
    
    departure_kep = cart2kep([departure_pos, departure_vel], constants.mu);
    
    departure_eq = kep2eq(departure_kep);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

% Equality constraints
Ceq = [];




GCeq = [];

GJ = [];

%% FABLE propagator 

% Final true longitude
L_f = arrival_eq(6);

% Initial conditions: six equinoctial elements, mass
E0 = [departure_eq; parameters.m];

if options.FABLE1
    % Variable to define position of variables in vector xopt;
    k = numel(x) - 4*arcs;
    % Vector of controls. It includes 4 columns.
    % 1st column: true longitude ON node
    % 2nd column: true longitude OFF node
    % 3rd column: azimuth angle
    % 4th column: elevation angle
    u(:,1) = x(3:4:end-k);
    u(:,2) = u(:,1) + x(4:4:end-k)';
    u(:,3) = x(1:4:end-k);
    u(:,4) = x(2:4:end-k)';
    
elseif options.FABLE2
    % Variable to define position of variables in vector xopt;
    k = numel(x) - 3*arcs;
    % Vector of controls. It includes 4 columns.
    % 1st column: acceleration
    % 2nd column: azimuth angle
    % 3rd column: elevation angle
    u(:,1) = x(1:3:end-k);
    u(:,2) = x(2:3:end-k);
    u(:,3) = x(3:3:end-k);
    
end

% Propagation
[ToF, x_final, J, output] = FABLE_propagator(L_f, E0, u, parameters, options, constants);

%% Objective and constraints

% Equality constraints
% %%%%%%%%%%%%%%% rendezvous %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ceq = [Ceq; x_final(1:5) - arrival_eq(1:5)];
% %%%%%%%%%%%%%%%%% fly-by %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% kep_final = eq2kep(x_final);
% cart_final = kep2cart(kep_final, constants.mu);
% 
% kep_arrival = eq2kep(arrival_eq);
% cart_arrival = kep2cart(kep_arrival, constants.mu);
% 
% Ceq = [Ceq; cart_final(1:3)'- cart_arrival(1:3)'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Objective
if ~options.J_DeltaV
    J_out = 0;
else
    J_out = J;
end

% Equality constraints.
k_ToF = 1e-3;
if options.const_ToF
    
    Ceq = [Ceq; (ToF - curr_tof)*k_ToF];
end

% Disequality non linear constraint - present only if the time has to be
% optimised
if options.const_ToF == 1
    C = [];
else
    % Final true longitude lower than arrival true longitude
    C = x(end - k  -1) + x(end - k ) - arrival_eq(6);
    % Sum of delta true longitude lower than difference between arrival and
    % departure true longitude
    C = [C; sum(x(4 : 4 : end - k )) - arrival_eq(6) + departure_eq(6)];
end

%% Gradient objective and constraints - check this, it might not be valid for FABLE2
if options.gradient == 1
    
    
    delta_FD = 1e-6;
    
    for ij = 1 : numel(x)

        x_var = x;
        x_var(ij) = x(ij) + delta_FD;
        
        % Equinoctial elements departure point
        if (numel(x) == 4*arcs && options.const_ToF == 1) || ...
            (numel(x) == 4*arcs + 1 && options.const_ToF == 0)
            % If departure state is given and declination and azimuth at departure
            % do not have to be optimised
            departure_eq = parameters.departure_eq;
        elseif (numel(x) == 4*arcs + 2 && options.const_ToF == 1) ||...
            (numel(x) == 4*arcs + 3 && options.const_ToF == 0) 
            % Departure position
            departure_pos = parameters.dep_body_r;
            % Velocity at departure is given by the sum of the velocity of the
            % central body and the velocity vinf
            v_inf = parameters.vinf0/ constants.DU * constants.TU * constants.sec_day;
            % Azimuth and declination
            alpha = x_var(end-k1-1);
            delta = x_var(end-k1);
%             if ij == 49
%                 keyboard
%             end
            departure_vel = parameters.dep_body_v + ...
                v_inf * [cos(alpha) * cos(delta), sin(alpha) * cos(delta), sin(delta)];
            
            departure_kep = cart2kep([departure_pos, departure_vel], constants.mu);
            
            departure_eq = kep2eq(departure_kep);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        
       
        % Equinoctial elements arrival point
        % 1st option: the time of flight is fixed, arrival point is defined
        if options.const_ToF == 1
            arrival_eq   = parameters.arrival_eq;
            L_f = arrival_eq(6);
        elseif options.const_ToF == 0 && ij == numel(x)
            % Compute arrival position based on time of fligth
            Apophis = inputs.Apophis;
            % 2nd option: recompute arrival position based on ToF
            % Adimensional semimajor axis
            kk = fix( ( (inputs.departure_date_MJD2000 + x_var(end)*constants.TU - Apophis.epochMJD2000 ) / constants.TU) ...
                /  Apophis.T );
            
            % Mean anomaly at arrival
            Apophis.M = mod(Apophis.M0 + ...
                Apophis.n * (x_var(end)*constants.TU + inputs.departure_date_MJD2000- Apophis.epochMJD2000)/constants.TU + ...
                -2 * kk * pi, 2*pi);
            
            % Eccentric anomaly at arrival
            Apophis.E = kepler_E(Apophis.e, Apophis.M);
            
            % True anomaly at arrival
            Apophis.cos_theta = ( cos(Apophis.E) - Apophis.e ) / ...
                (1 - Apophis.e * cos(Apophis.E));
            
            Apophis.sin_theta = ( sin(Apophis.E)  * sqrt(1 - Apophis.e^2) ) / ...
                (1 - Apophis.e * cos(Apophis.E));
            
            Apophis.theta = atan2(Apophis.sin_theta,Apophis.cos_theta);
            
            % Orbital element of apophis at arrival
            arrival_kep = [Apophis.a, Apophis.e, Apophis.i, Apophis.RAAN, Apophis.omega, ...
                Apophis.theta];
            
            arrival_eq = kep2eq(arrival_kep);
            
            % Increase the arrival true longitude by 2*pi times the number of
            % revolutions
            arrival_eq(6) = arrival_eq(6) + 2 * pi * (options.n_rev);
            L_f = arrival_eq(6);
        end


        % Initial conditions: six equinoctial elements, mass
        E0 = [departure_eq; parameters.m];
        
        % Vector of controls. It includes 4 columns.
        % 1st column: true longitude ON node
        % 2nd column: true longitude OFF node
        % 3rd column: azimuth angle
        % 4th column: elevation angle
        u_var(:,1) = x_var(3:4:end-k);
        u_var(:,2) = u_var(:,1) + x_var(4:4:end-k)';
        u_var(:,3) = x_var(1:4:end-k);
        u_var(:,4) = x_var(2:4:end-k)';
        
        % Propagation
        [ToF_var, x_final_var, J_var] = FABLE_propagator(L_f, E0, u_var, parameters, options, constants);
        
        % Gradient equality constraints
        %%%%%%%%%%%%%%%%%%%%%% rendezvous %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        GCeq(1:5,ij) = (x_final_var(1:5) - x_final(1:5) ) / delta_FD;
        if options.const_ToF == 1
            GCeq(6, ij)  = (ToF_var - ToF)*k_ToF / delta_FD;
        end
        %%%%%%%%%%%%%%%%%%%%%% flyby %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         kep_final_var = eq2kep(x_final_var);
%         cart_final_var = kep2cart(kep_final_var, constants.mu);
%         
% %         kep_arrival = eq2kep(arrival_eq);
% %         cart_arrival = kep2cart(kep_arrival, constants.mu);
%         
%         GCeq(1:3,ij) = (cart_final_var(1:3) - cart_final(1:3) ) / delta_FD;
%         if options.const_ToF == 1
%             GCeq(4, ij)  = (ToF_var - ToF)*k_ToF / delta_FD;
%         end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        % Gradient objective function
        GJ(ij)  = (J_var - J ) / delta_FD;
        
        if options.const_ToF == 0
            % Gradient inequality constraints
            GC(1,ij) = [x_var(end - k  -1) + x_var(end - k ) - arrival_eq(6) - C(1)] / delta_FD;
            GC(2,ij) = [ sum(x_var(4 : 4 : end - k )) - arrival_eq(6) + departure_eq(6) - C(2)] / delta_FD;
        else
            GC = [];
        end
        
        
    end
    
    
else
    GC = [];
end



%% Other outputs

Equin_thrust = output.Equin_thrust;
Equin_coast = output.Equin_coast;

Equin_sorted = output.Equin_sorted;
Equin_sorted_ct = output.Equin_sorted_ct;
time_sorted = output.time_sorted;

%% Plot

if parameters.plot_flag
    
    % Equinoctial elements departure point
    if (numel(x) == 4*arcs && options.const_ToF == 1) || ...
            (numel(x) == 4*arcs + 1 && options.const_ToF == 0)
        % If departure state is given and declination and azimuth at departure
        % do not have to be optimised
        departure_eq = parameters.departure_eq;
    elseif (numel(x) == 4*arcs + 2 && options.const_ToF == 1) ||...
            (numel(x) == 4*arcs + 3 && options.const_ToF == 0)
        % Departure position
        departure_pos = parameters.dep_body_r;
        % Velocity at departure is given by the sum of the velocity of the
        % central body and the velocity vinf
        v_inf = parameters.vinf0/ constants.DU * constants.TU * constants.sec_day;
        % Azimuth and declination
        alpha = x(end-k1-1);
        delta = x(end-k1);
        %             if ij == 49
        %                 keyboard
        %             end
        departure_vel = parameters.dep_body_v + ...
            v_inf * [cos(alpha) * cos(delta), sin(alpha) * cos(delta), sin(delta)];
        
        departure_kep = cart2kep([departure_pos, departure_vel], constants.mu);
        
        departure_eq = kep2eq(departure_kep);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    
    % Equinoctial elements arrival point
    % 1st option: the time of flight is fixed, arrival point is defined
    if options.const_ToF == 1
        arrival_eq   = parameters.arrival_eq;
        L_f = arrival_eq(6);
    elseif options.const_ToF == 0 && ij == numel(x)
        % Compute arrival position based on time of fligth
        Apophis = inputs.Apophis;
        % 2nd option: recompute arrival position based on ToF
        % Adimensional semimajor axis
        kk = fix( ( (inputs.departure_date_MJD2000 + x(end)*constants.TU - Apophis.epochMJD2000 ) / constants.TU) ...
            /  Apophis.T );
        
        % Mean anomaly at arrival
        Apophis.M = mod(Apophis.M0 + ...
            Apophis.n * (x(end)*constants.TU + inputs.departure_date_MJD2000- Apophis.epochMJD2000)/constants.TU + ...
            -2 * kk * pi, 2*pi);
        
        % Eccentric anomaly at arrival
        Apophis.E = kepler_E(Apophis.e, Apophis.M);
        
        % True anomaly at arrival
        Apophis.cos_theta = ( cos(Apophis.E) - Apophis.e ) / ...
            (1 - Apophis.e * cos(Apophis.E));
        
        Apophis.sin_theta = ( sin(Apophis.E)  * sqrt(1 - Apophis.e^2) ) / ...
            (1 - Apophis.e * cos(Apophis.E));
        
        Apophis.theta = atan2(Apophis.sin_theta,Apophis.cos_theta);
        
        % Orbital element of apophis at arrival
        arrival_kep = [Apophis.a, Apophis.e, Apophis.i, Apophis.RAAN, Apophis.omega, ...
            Apophis.theta];
        
        arrival_eq = kep2eq(arrival_kep);
        
        % Increase the arrival true longitude by 2*pi times the number of
        % revolutions
        arrival_eq(6) = arrival_eq(6) + 2 * pi * (options.n_rev);
        L_f = arrival_eq(6);
    end
    
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
    
end

