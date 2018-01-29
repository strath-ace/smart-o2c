function [DeltaV, time, xopt, Equin_sorted, Equin_sorted_ct, time_sorted, ...
    LowThrustTime, parameters, GCeq, GJ] = FABLE ...
    (departure_kep, arrival_kep, ToF, engine, spacecraft, options, constants, inputs)

% =========================================================================
% FABLE: Fast Analytic Boundary Value Estimator.
% 
% Direct optimisation of low-thrust transfers using Multiple Shooting (MS)
% or Single Shooting (SS) and Matlab fmincon (either SQP or IPM).
% =========================================================================
% Input: departure_kep -> keplerian elements of the departure point
%                         (if the v_inf vector at departure is optimised,
%                         this input is discarded and recomputed)
%        arrival_kep   -> keplerian elements of the arrival point
%        ToF           -> time of flight for the transfer
%        engine        -> structure with information about the low-thrust
%                         engine (thrust, acceleration, specific impulse)
%        spacecraft    -> structure with information about the spacecrat
%        options       -> options for the transfer (number of transfer
%                         arcs, IG provided by the user)
%        constants     -> structure with constants for the problem
%        inputs        -> (inputs to the problem)
% 
% Output: DeltaV          -> Delta V for the transfer
%         time            -> time of the transfer
%         xopt            -> optimal solution vector
%         Equin_sorted    -> cells of data. Each cell contains equinoctial
%                            elements and mass variation of each
%                            thrust/coast arc
%         Equin_sorted_ct -> vector of 0s and 1s. 1 when the engine is on,
%                            0 when the engine is off
%         time_sorted     -> cells of data for the time. Same structure as
%                            Equin_sorted
%         LowThrustTime   -> time with engine on
%         parameters      -> 
%         GCeq            -> matrix of gradients of the constraints
%         GJ              -> gradients of the objectives
% =========================================================================
% References: 
% 1. F. Zuiani, 
% 2. 
% =========================================================================
% Marilena Di Carlo, marilena.di-carlo@strath.ac.uk

%% Inputs

% Number of thrust arcs
arcs = options.transfer_arcs;

engine.acc = engine.T_adim/ spacecraft.m;

% Acceleration
eps_max = engine.acc;

% Isp
parameters.Isp_adim = engine.Isp_adim;

% Low-thrust acceleration [DU/TU^2]
% eps_max = eps_max * 10^(-3) * (((constants.TU*constants.sec_day)^2) / constants.DU);


%% Departure and arrival equinoctial elements

% If the state at derarture is to be optimised, recompute departure_kep,
% otherwise use the one defined in the input
if isfield(inputs,'v_inf')
    
    % Position and velocity of the departure body at departure
    parameters.dep_body_r = inputs.Earth_r;
    parameters.dep_body_v = inputs.Earth_v;
    
    % Magnitude of the vinf vector
    vinf0 = inputs.v_inf;
    parameters.vinf0 = vinf0;
    
    % Position of the departure
    departure_pos = parameters.dep_body_r;
    
    % Velocity at departure: velocity of the departure body + v_inf
    % velocity (TO DO: v_inf velocity should be defined in the geocentric
    % equatorial plane. Here it is defined in a geocentric ecliptic plane.
    % alpha and delta can however be rotated to obtain angles in geocentric
    % equatorial plane)
    departure_vel = parameters.dep_body_v + ...
        vinf0/ constants.DU * constants.TU * constants.sec_day *...
        [cos(inputs.alpha0) * cos(inputs.delta0), ...
        sin(inputs.alpha0) * cos(inputs.delta0), ...
        sin(inputs.delta0)];
    
    % Keplerian elements at departure
    departure_kep = cart2kep([departure_pos, departure_vel], constants.mu);
    
    % Equinoctial elements at departure
    departure_eq = kep2eq(departure_kep);
    
else
    % Equinoctial elements at departure
    departure_eq = kep2eq(departure_kep);
    
end


% Equinoctial elements of the ARRIVAL point
arrival_eq = kep2eq(arrival_kep);



%% True longitude variation

% Only counterclockwise motion of the spacecraft
if arrival_eq(6) < departure_eq(6)
    
    arrival_eq(6) = arrival_eq(6) + 2*pi;
    
end

% Increase the arrival true longitude by 2*pi times the number of
% revolutions
arrival_eq(6) = arrival_eq(6) + 2 * pi * (options.n_rev);


%% Initial guess vector

% If the initial guess is not given as input
if options.IG_flag == 0
    
    
    % Initial guess for multiple shooting
    if strcmp(options.MS_SS,'MS')
        
        if options.FABLE1
            n_eq_nodes = 2*arcs;
            % Variation of true longitude from departure to arrival point
            dL = linspace(departure_eq(6), arrival_eq(6), 2*arcs+2);
        elseif options.FABLE2
            n_eq_nodes = arcs-1;
        end
        
        % The vector x is defined by 14 elements per arc.
        % x(1)     - azimuth angle for thrust vector
        % x(2)      - elevation angle
        % x(3:8)   - orbital elements of the initial point of the 1st thrust arc
        %           (a, P1, P2, Q1, Q2, L)
        % x(9:14)  - orbital elements of the intial poinf of the coast arc
        %           (a, P1 ,P2, Q1, Q2, Delta L)
        % where Delta L is the difference in L wrt to L (x(8))
        
        

        
        % Nodes of the shooting method
        if strcmp(options.first_guess, 'departure')
            equin_nodes(1,:) = linspace(departure_eq(1), arrival_eq(1), n_eq_nodes);
            equin_nodes(2,:) = departure_eq(2) * ones(1,n_eq_nodes);
            equin_nodes(3,:) = departure_eq(3) * ones(1,n_eq_nodes);
            equin_nodes(4,:) = departure_eq(4) * ones(1,n_eq_nodes);
            equin_nodes(5,:) = departure_eq(5) * ones(1,n_eq_nodes);
            
        elseif strcmp(options.first_guess, 'arrival')
            
            equin_nodes(1,:) = linspace(departure_eq(1), arrival_eq(1), n_eq_nodes);
            equin_nodes(2,:) = arrival_eq(2) * ones(1,n_eq_nodes);
            equin_nodes(3,:) = arrival_eq(3) * ones(1,n_eq_nodes);
            equin_nodes(4,:) = arrival_eq(4) * ones(1,n_eq_nodes);
            equin_nodes(5,:) = arrival_eq(5) * ones(1,n_eq_nodes);
            
        elseif strcmp(options.first_guess, 'linspace')
            %
            equin_nodes(1,:) = linspace(departure_eq(1), arrival_eq(1), n_eq_nodes);
            equin_nodes(2,:) = linspace(departure_eq(2), arrival_eq(2), n_eq_nodes);
            equin_nodes(3,:) = linspace(departure_eq(3), arrival_eq(3), n_eq_nodes);
            equin_nodes(4,:) = linspace(departure_eq(4), arrival_eq(4), n_eq_nodes);
            equin_nodes(5,:) = linspace(departure_eq(5), arrival_eq(5), n_eq_nodes);
            
            
        end
        
        % -----------Initial guess FABLE 1 --------------------------------
        if options.FABLE1
            for i = 1 : size(equin_nodes,2)
                equin_nodes(6,i) = dL(i+1);
            end
            
            x00 = zeros(arcs, 14);
            
            % Azimuth angle thrust vector
            if arrival_eq(1) >= departure_eq(1)
                x00(:,1) = 90*(pi/180);
            else
                x00(:,1) = -90*(pi/180);
            end
            
            % Elevation angle
            x00(:,2) = 0*(pi/180);
            
            for i = 1 :  arcs
                
                x00(i,3:8)  = [equin_nodes(1:5,1+(i-1)*2)' equin_nodes(6,1+(i-1)*2)];
                x00(i,9:14)  = [equin_nodes(1:5,2*i)' equin_nodes(6,2*i)-equin_nodes(6,1+(i-1)*2)];
                
            end
            
            x0 = reshape(x00', 1, 14*arcs);
            
            % -----------End initial guess FABLE 1 ----------------------------
        elseif options.FABLE2
            % -----------Initial guess FABLE 2 ----------------------------
            x00 = zeros(arcs-1, 8);
            
            % Magnitude acceleration
            x00(:,1) = options.T0*(eps_max - 0.01*eps_max);
            
            % Azimuth angle
            if arrival_eq(1) >= departure_eq(1)
                x00(:,2) = 90*(pi/180);
            else
                x00(:,2) = -90*(pi/180);
            end
            
            % Elevation angle
            x00(:,3) = 0*(pi/180);
            
            for i = 1 :  arcs-1
                
                x00(i,4:8)  = equin_nodes(1:5,i)';
                
            end
            
            % Reshape vector
            x0 = reshape(x00', 1, 8*(arcs-1));
            
            x0 = [options.T0*(eps_max- 0.01*eps_max), x00(1,2), 0, x0];
        end
        
        
 
        
        % If the departure vinf vector is to be optimised in azimuth and
        % declination, add variables
        if isfield(inputs,'v_inf')
            x0 = [x0,  inputs.alpha0, inputs.delta0];
        end
        
        
        % Initial guess for single shooting
    elseif strcmp(options.MS_SS,'SS')
        
        if options.FABLE1
            % Variation of true longitude from departure to arrival point
            dL = linspace(departure_eq(6), arrival_eq(6), 2*arcs+2);
            
            for i = 1 : 2*arcs
                equin_nodes(6,i) = dL(i+1);
            end
            
            x00 = zeros(arcs, 4);
            
            for i = 1 :  arcs
                
                x00(i,3)  = equin_nodes(6,1+(i-1)*2);
                x00(i,4)  = equin_nodes(6,2*i)-equin_nodes(6,1+(i-1)*2);
                
            end
            
            % Azimuth angle thrust vector
            if arrival_eq(1) >= departure_eq(1)
                x00(:,1) = 90*(pi/180);
            else
                x00(:,1) = -90*(pi/180);
            end
            
            % Elevation angle
            x00(:,2) = 0*(pi/180);
            
            x0 = reshape(x00', 1, 4*arcs);
            
            
        elseif options.FABLE2
            
            x00 = zeros(arcs, 3);
            
            % Magnitude acceleration
            x00(:,1) = options.T0*(eps_max - 0.01*eps_max) / spacecraft.m;
            
            % Azimuth angle thrust vector
            if arrival_eq(1) >= departure_eq(1)
                x00(:,2) = 90*(pi/180);
            else
                x00(:,2) = -90*(pi/180);
            end
            
            % Elevation angle
            x00(:,3) = 0*(pi/180);
            
            x0 = reshape(x00', 1, 3*arcs);
            
        end
        
        % If the departure vinf vector is to be optimised in azimuth and
        % declination, add variables
        if isfield(inputs,'v_inf')
            x0 = [x0,  inputs.alpha0, inputs.delta0];
        end
        
        % If the time of flight has to be optimised, add variable
        if options.const_ToF == 0
            x0 = [x0, inputs.ToF];
        end
        
    end
    
else
    
    x0 = inputs.IG;
end

%% Lower and upper boundaries

min_semimajor = min(departure_eq(1), arrival_eq(1));
max_semimajor = max(departure_eq(1), arrival_eq(1));

min_semimajor = min_semimajor - 0.01*min_semimajor;
max_semimajor = max_semimajor + 0.01*max_semimajor;
% min_semimajor = 0.8 * constants.AU / constants.DU;
% max_semimajor = 1.2* constants.AU / constants.DU;


max_incl = max(departure_kep(3), arrival_kep(3));
max_incl = max_incl + 0.1 * max_incl;

Q1max = tan(max_incl / 2);


% Boundaries for multiple shooting
if strcmp(options.MS_SS,'MS')
    
    if options.FABLE1
        LB = reshape([-180* (pi/180);  -90*pi/180;  ...
            min_semimajor - 0.1 * min_semimajor;  -1; -1;   -Q1max; -Q1max; departure_eq(6);  ...
            min_semimajor - 0.1 * min_semimajor;  -1; -1;   -Q1max; -Q1max; 0] * ...
            ones(1,arcs) , 1, arcs+13*arcs);
        
        UB = reshape([ 180* (pi/180);   90*pi/180;  ...
            max_semimajor + 0.1 * max_semimajor;   1;  1;   Q1max; Q1max; arrival_eq(6); ...
            max_semimajor + 0.1 * max_semimajor;   1;  1;   Q1max; Q1max; ...
            arrival_eq(6) - departure_eq(6)]   ...
            * ones(1,arcs) , 1, arcs+13*arcs);
        
    elseif options.FABLE2
        
        LB = reshape([0; -180* (pi/180);  -90*pi/180;  ...
            min_semimajor - 0.1 * min_semimajor;  -1; -1;   -Q1max; -Q1max] *...
            ones(1,arcs-1) , 1, 8*(arcs-1));
        UB = reshape([eps_max;  180* (pi/180);   90*pi/180;  ...
            max_semimajor + 0.1 * max_semimajor;   1;  1;   Q1max; Q1max]   ...
            * ones(1,arcs-1) , 1, 8*(arcs-1));
        
        LB = [0, -180* (pi/180),  -90*pi/180, LB];
        UB = [eps_max,  180* (pi/180),   90*pi/180, UB];
        
    end
    
    if isfield(inputs,'v_inf')
        LB = [LB, 0, -pi/2];
        UB = [UB,  2*pi, 0];
    end
    
    % Boundaries for single shooting
elseif  strcmp(options.MS_SS,'SS')
    
    if options.FABLE1
        
        LB = reshape([-180* (pi/180);  -90*pi/180;  ...
            departure_eq(6);  ...
            0] * ...
            ones(1,arcs) , 1, 4*arcs);
        
        UB = reshape([ 180* (pi/180);   90*pi/180;  ...
            arrival_eq(6); ...
            Inf]   ...
            * ones(1,arcs) , 1, 4*arcs);
        
    elseif options.FABLE2
        LB = reshape([0; -180* (pi/180);  -90*pi/180] * ...
            ones(1,arcs) , 1, 3*arcs);
        
        UB = reshape([eps_max;  180* (pi/180);   90*pi/180]   ...
            * ones(1,arcs) , 1, 3*arcs);
    end
    
    if isfield(inputs,'v_inf')
        LB = [LB, 0, -pi/2];
        UB = [UB,  2*pi, 0];
    end
    
    if options.const_ToF == 0
        LB = [LB, inputs.ToF_LB];
        UB = [UB, inputs.ToF_UB];
    end
    
end

%% Linear constraints

% INEQUALITIES CONSTRAINTS ON TRUE LONGITUDE VALUES
% Inequalities conditions:
% departure_eq(6) < L1 < L2 < L3 and so on until Ln < arrival_eq(6)
% Inequalities constraints are expressed through the relationship Ax<b


% Multiple shooting
if strcmp(options.MS_SS,'MS')
    
    if isfield(inputs,'v_inf')
        A = zeros(2*arcs, 14*arcs+2);
    else
        A = zeros(2*arcs, 14*arcs);
    end
    b = zeros(2*arcs,1);
    for i = 1 : arcs-1
        A(i ,8 + 14 * (i-1)) = 1;
        A(i, 8+14+ 14 * (i-1)) = -1;
    end
    for i = 1 : arcs-1
        A(arcs-1+i,8 + 14 * (i-1)) = 1;
        A(arcs-1+i,14 + 14 * (i-1)) = 1;
        A(arcs-1+i,8 +14+ 14 * (i-1)) = -1;
    end
    A(end-1,8 + 14 * (i)) = 1;
    A(end-1,14 + 14 * (i)) = 1;
    b(end-1,1) = arrival_eq(6);
    A(end,14:14:end) = 1;
%     b(end) = 2*pi*options.n_rev;
    b(end) = arrival_eq(6) - departure_eq(6);
    
    Aeq = [];
    beq = [];
    
    
    % Single shooting
elseif strcmp(options.MS_SS,'SS')
    
    if isfield(inputs,'v_inf') && options.const_ToF == 1
        A = zeros(2*arcs, 4*arcs+2);
    elseif isfield(inputs,'v_inf') && options.const_ToF == 0
        A = zeros(2*arcs, 4*arcs+3);
    elseif ~isfield(inputs,'v_inf') && options.const_ToF == 0
        A = zeros(2*arcs, 4*arcs+1);
    else
        A = zeros(2*arcs, 4*arcs);
    end
    
    b = zeros(2*arcs,1);
    
    for i = 1 : arcs-1
        A(i ,3 + 4 * (i-1)) = 1;
        A(i, 7+ 4 * (i-1)) = -1;
    end
    
    for i = 1 : arcs-1
        A(arcs-1+i,3 + 4 * (i-1)) = 1;
        A(arcs-1+i,4 + 4 * (i-1)) = 1;
        A(arcs-1+i,7+ 4 * (i-1)) = -1;
    end
    
    if options.const_ToF == 1
        A(end-1,3 + 4 * (i)) = 1;
        A(end-1,4 + 4 * (i)) = 1;
        b(end-1,1) = arrival_eq(6);
        A(end,4:4:end) = 1;
        b(end) = 2*pi*options.n_rev;
        b(end) = arrival_eq(6) - departure_eq(6);
    end
    
    Aeq = [];
    beq = [];
    
end



if options.FABLE2
    A = [];
    b = [];
end

%% Optimization

parameters.arcs         = arcs;
parameters.eps_max      = eps_max;
parameters.departure_eq = departure_eq;
parameters.arrival_eq   = arrival_eq;
parameters.curr_tof     = ToF;
parameters.plot_flag    = 0;
parameters.n            = 2;
parameters.m            = spacecraft.m;
parameters.T            = eps_max * spacecraft.m;
parameters.LB           = LB;
parameters.UB           = UB;

% Check if zero thrust solution is the solution
if strcmp(options.MS_SS,'MS') 
    [~,~,Ceq_x0] = FABLE_transcription_MS(x0, parameters, options, constants);
elseif strcmp(options.MS_SS,'SS')
    [~,~,Ceq_x0] = FABLE_transcription_SS(x0, parameters, options, constants, inputs);
end
x0_C_viol = max(abs(Ceq_x0));

if x0_C_viol == 0
    DeltaV = 0;
    time = 0;
    xopt = [];
    Equin_sorted = [];
    Equin_sorted_ct=[];
    
else
    
    options.options_fmincon = optimset('Display',options.options_fmincon_temp.Display,...
        'MaxFunEvals',options.options_fmincon_temp.MaxFunEvals,...
        'LargeScale',options.options_fmincon_temp.LargeScale,...
        'Algorithm',options.options_fmincon_temp.Algorithm,...
        'TolCon', options.options_fmincon_temp.TolCon ,...
        'TolX',options.options_fmincon_temp.TolX, ...
        'GradConstr',options.options_fmincon_temp.GradConstr,...
        'GradObj',options.options_fmincon_temp.GradObj,...
        'DerivativeCheck',options.options_fmincon_temp.DerivativeCheck, ...
        'MaxIter',options.options_fmincon_temp.MaxIter);
    
    % Optimisation
    [xopt,fopt,iout,output] = FABLE_runobjconstr(x0, parameters, options, ...
        constants, inputs, A, b, Aeq, beq,  LB, UB);
    
    % Change values of the parameters to allow for the plot
    parameters.plot_flag = options.plot_flag;
    parameters.n         = 100;
    
    if strcmp(options.MS_SS,'MS') 
        [J_opt,C_opt,Ceq_opt2,GJ,GC,GCeq,ToF,Equin_sorted,Equin_sorted_ct,time_sorted,J_DV] = ...
            FABLE_transcription_MS(xopt, parameters, options, constants);
    elseif strcmp(options.MS_SS,'SS')
        [J_opt,C_opt,Ceq_opt2,GJ, GC, GCeq,ToF,Equin_sorted,Equin_sorted_ct,time_sorted,J_DV] = ...
            FABLE_transcription_SS(xopt, parameters, options, constants, inputs);
    end
    
    time = ToF * constants.TU;
    
    
    
    if (iout == 1  ||  iout == 2)
        
        if options.fun == 0
            DeltaV = J_DV * eps_max / ((constants.sec_day*constants.TU)/constants.DU);
        elseif options.fun == 1
            % deltaV in km/s (AU is in km)
            DeltaV = J_DV  /  ((constants.sec_day*constants.TU)/constants.DU);
        end
        
        % Time spent using low thrust propulsion
        LowThrustTime = 0;
        for i = 1 : numel(Equin_sorted_ct)
            if ~isempty(time_sorted{i})
                LowThrustTime = LowThrustTime + Equin_sorted_ct(i) .* ...
                ( time_sorted{i}(end) - time_sorted{i}(1) );
            end
        end
        
        
        
    else
        DeltaV = NaN;
        LowThrustTime = NaN;
    end
    
end



end