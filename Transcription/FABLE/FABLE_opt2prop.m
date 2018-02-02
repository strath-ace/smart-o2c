function [L_f, E0, u, parameters_prop, arrival_eq] = FABLE_opt2prop(xopt, departure_eq, ...
    arrival_eq, engine, spacecraft, options, inputs, constants)

%% Final true longitude

if options.const_ToF == 1
    L_f = arrival_eq(6) + 2 * pi * (options.n_rev);
else
     Apophis = inputs.Apophis;
    % 2nd option: recompute arrival position based on ToF
    % Adimensional semimajor axis    
    kk = fix( ( (inputs.departure_date_MJD2000 + xopt(end)*constants.TU - Apophis.epochMJD2000 ) / constants.TU) ...
        /  Apophis.T );
    
    % Mean anomaly at arrival
    Apophis.M = mod(Apophis.M0 + ...
        Apophis.n * (xopt(end)*constants.TU + inputs.departure_date_MJD2000- Apophis.epochMJD2000)/constants.TU + ...
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
    
    L_f = arrival_eq(6) ;
end


%% Initial state

% Variable to define position of variables in vector xopt;
k1 = options.const_ToF ==  0;

if strcmp(options.MS_SS,'MS')
    
    % Equinoctial elements departure point
    if options.FABLE1 && numel(xopt) == 14*options.transfer_arcs || ...
            options.FABLE2 && numel(xopt) == 8*(options.transfer_arcs-1)+3
        % If departure state is given and declination and azimuth at departure
        % do not have to be optimised
        departure_eq = departure_eq;
    elseif options.FABLE1 && numel(xopt) == 14*options.transfer_arcs + 2 || ...
            options.FABLE2 && numel(xopt) == 8*(options.transfer_arcs-1)+3 + 2
        % Departure position
        departure_pos = inputs.Earth_r;
        % Velocity at departure is given by the sum of the velocity of the
        % central body and the velocity vinf
        v_inf = inputs.v_inf/ constants.DU * constants.TU * constants.sec_day;
        % Azimuth and declination
        alpha = xopt(end-1);
        delta = xopt(end);
        
        departure_vel = inputs.Earth_v + ...
            v_inf * [cos(alpha) * cos(delta), sin(alpha) * cos(delta), sin(delta)];
        
        departure_kep = cart2kep([departure_pos, departure_vel], constants.mu);
        
        departure_eq = kep2eq(departure_kep);
    end
    
    
elseif strcmp(options.MS_SS,'SS')
    
    % Equinoctial elements departure point
    if (options.FABLE1 && numel(xopt) == 4*options.transfer_arcs  && options.const_ToF == 1) || ...
       (options.FABLE1 && numel(xopt) == 4*options.transfer_arcs+1  && options.const_ToF == 0) || ...
        options.FABLE2 && numel(xopt) == 3*(options.transfer_arcs-1)+3
        % If departure state is given and declination and azimuth at departure
        % do not have to be optimised
        departure_eq = departure_eq;
    elseif (options.FABLE1 && numel(xopt) == 4*options.transfer_arcs + 2 && options.const_ToF == 1) || ...
           (options.FABLE1 && numel(xopt) == 4*options.transfer_arcs + 3 && options.const_ToF == 0) || ...
            options.FABLE2 && numel(xopt) == 3*(options.transfer_arcs-1)+3 + 2
        % Departure position
        departure_pos = inputs.Earth_r;
        % Velocity at departure is given by the sum of the velocity of the
        % central body and the velocity vinf
        v_inf = inputs.v_inf/ constants.DU * constants.TU * constants.sec_day;
        % Azimuth and declination
        alpha = xopt(end-k1-1);
        delta = xopt(end-k1);
        
        departure_vel = inputs.Earth_v + ...
            v_inf * [cos(alpha) * cos(delta), sin(alpha) * cos(delta), sin(delta)];
        
        departure_kep = cart2kep([departure_pos, departure_vel], constants.mu);
        
        departure_eq = kep2eq(departure_kep);
    end
    
end


    
% Initial conditions: six equinoctial elements, mass
E0 = [departure_eq; spacecraft.m];

%% Controls

% Control matrix
if strcmp(options.MS_SS,'MS')
    if options.FABLE1
        % Variable to define position of variables in vector xopt;
        k = numel(xopt) - 14*options.transfer_arcs;
        % Vector of controls. It includes 4 columns.
        % 1st column: true longitude ON node
        % 2nd column: true longitude OFF node
        % 3rd column: azimuth angle
        % 4th column: elevation angle
        u(:,1) = xopt(8:14:end-k);
        u(:,2) = u(:,1) + xopt(14:14:end-k)';
        u(:,3) = xopt(1:14:end-k)';
        u(:,4) = xopt(2:14:end-k)';
    elseif options.FABLE2
        % Variable to define position of variables in vector xopt;
        k = numel(xopt) - 8*(options.transfer_arcs-1)-3;
        % Vector of controls. It includes 3 columns.
        % 1st column: acceleration
        % 2nd column: azimtuh angle
        % 3rd column: elevation angle
        u(1,1:3) = xopt(1:3);
        u = [u; xopt(4 : 8 : end - k)' xopt(5 : 8 : end - k)' ...
            xopt(6 : 8 : end - k)'];
        %         u(2:end,1) = xopt(4 : 8 : end - k)';
        %         u(2:end,2) = xopt(5 : 8 : end - k)';
        %         u(3:end,2) = xopt(6 : 8 : end - k)';
    end
    
elseif strcmp(options.MS_SS,'SS')
    if options.FABLE1
        % Variable to define position of variables in vector xopt;
        k = numel(xopt) - 4*options.transfer_arcs;
        
        % Vector of controls. It includes 4 columns.
        % 1st column: true longitude ON node
        % 2nd column: true longitude OFF node
        % 3rd column: azimuth angle
        % 4th column: elevation angle
        u(:,1) = xopt(3:4:end-k);
        u(:,2) = u(:,1) + xopt(4:4:end-k)';
        u(:,3) = xopt(1:4:end-k)';
        u(:,4) = xopt(2:4:end-k)';
        
    elseif options.FABLE2
        % Variable to define position of variables in vector xopt;
        k = numel(xopt) - 3*options.transfer_arcs;
        % Vector of controls. It includes 3 columns.
        % 1st column: acceleration
        % 2nd column: azimtuh angle
        % 3rd column: elevation angle
        u(:,1) = xopt(1 : 3 : end - k);
        u(:,2) = xopt(2 : 3 : end-k);
        u(:,3) = xopt(3:3:end-k);
        
    end
end

%% Other propagation parameters

parameters_prop.Isp_adim = engine.Isp_adim;
parameters_prop.arcs   = options.transfer_arcs;
parameters_prop.m      = spacecraft.m;
parameters_prop.T_adim = engine.T_adim;
parameters_prop.n      = 100;

