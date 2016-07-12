classdef KeplerianElements
     
    properties
        
        a    % Semimajor axis [AU]     
        e    % Eccentricity
        i    % Inclination [deg]     
        OM   % Asc. Node/raan [deg]
        W    % Arg. Perigee [deg]     
        M0   % Mean anomoly, M at time given t0 [deg]
        t0   % Time at which Mo is given [MJD2000]   
        
    end
    
    
    methods (Access = public)
                       
        % Default KeplerianElements Constructor
        %  
        % Params: 
        %  a   Semimajor axis [AU]     
        %  e   Eccentricity 
        %  i   Inclination [deg]     
        %  OM  Asc. Node/raan [deg]
        %  W   Arg. Perigee [deg]     
        %  M0  Mean anomoly, M at time given t0 [deg]
        %  t0  Time at which Mo is given [MJD2000]  
        function kep_obj = KeplerianElements(a,e,i,om,w,M0,t0) 
           kep_obj.a  = a; 
           kep_obj.e  = e;
           kep_obj.i  = i;
           kep_obj.OM = om;
           kep_obj.W  = w;
           kep_obj.M0 = M0;
           kep_obj.t0 = t0;

        end
           
        
        % Get the Orbital Elements Array including True Anomaly in DEGREES
        %
        %  elementsarray = [a e i Om om theta]    
        %
        %  Where:
        %     a      Semimajor axis [AU]     
        %     e      Eccentricity 
        %     i      Inclination [deg]     
        %     OM     Asc. Node/raan [deg]
        %     om     Arg. Perigee [deg]     
        %     theta  True anomaly [deg]  
        %
        function elementsarray = getKepArray_AU_DEG(obj)
            
            % Get the True Anomaly
            f = getTrueAnomaly(obj);
                       
            % Create the Orbital Elements Arrays [a e i Om om theta]            
            elementsarray = [obj.a,  ... % 
                             obj.e,  ... % 
                             obj.i,  ... % 
                             obj.OM, ... % 
                             obj.W,  ... %  
                             f];     ... % True Anomaly 
            
        end
                   
        
        % Get the Orbital Elements Array including True Anomaly in RADIANS
        %
        %  elementsarray = [a e i Om om theta]    
        %
        %  Where:
        %     a      Semimajor axis [AU]     
        %     e      Eccentricity 
        %     i      Inclination [rad]     
        %     OM     Asc. Node/raan [rad]
        %     om     Arg. Perigee [rad]     
        %     theta  True anomaly [rad]
        %
        function elementsarray = getKepArray_AU_RAD(obj)
                        
            %  Get Keplerian Elements Array Including True Anomaly in DEG
            temp = obj.getKepArray_AU_DEG();
                        
            % Create the Orbital Elements Arrays [a e i Om om theta] [AU, RAD]                   
            elementsarray = [temp(1),          ... %
                             temp(2),          ... %
                             temp(3) * pi/180, ... %
                             temp(4) * pi/180, ... %
                             temp(5) * pi/180, ... %
                             temp(6) * pi/180];    %
                         
        end
            
        
           % Get the Orbital Elements Array including True Anomaly in DEGREES
        %
        %  elementsarray = [a e i Om om theta]    
        %
        %  Where:
        %     a      Semimajor axis [km]     
        %     e      Eccentricity 
        %     i      Inclination [deg]     
        %     OM     Asc. Node/raan [deg]
        %     om     Arg. Perigee [deg]     
        %     theta  True anomaly [deg]  
        %
        function elementsarray = getKepArray_KM_DEG(obj)
                       
            % Convertsion factor
            global Astroconstants
            au2km = AstroConstants.Astronomical_Unit;
            
            % Get the True Anomaly
            f = getTrueAnomaly(obj);
            
            % Create the Orbital Elements Arrays [a e i Om om theta]            
            elementsarray = [obj.a * au2km, ... % 
                             obj.e,         ... %
                             obj.i,         ... %
                             obj.OM,        ... %
                             obj.W,         ... %
                             f];
            
        end
                   
        
        % Get the Orbital Elements Array including True Anomaly in RADIANS
        %
        %  elementsarray = [a e i Om om theta]    
        %
        %  Where:
        %     a      Semimajor axis [km]     
        %     e      Eccentricity 
        %     i      Inclination [rad]     
        %     OM     Asc. Node/raan [rad]
        %     om     Arg. Perigee [rad]     
        %     theta  True anomaly [rad]
        %
        function elementsarray = getKepArray_KM_RAD(obj)
                        
            % Convertsion factor
            au2km = AstroConstants.Astronomical_Unit;
            
            %  Get Keplerian Elements Array Including True Anomaly in DEG
            temp = obj.getKepArray_AU_DEG();                        
            
            % Create the Orbital Elements Arrays [a e i Om om theta] [AU, RAD]                   
            elementsarray = [temp(1) * au2km,  ... % Semimajor [km]
                             temp(2),          ... % Eccentricity 
                             temp(3) * pi/180, ... % Inclination [rad]
                             temp(4) * pi/180, ... % Asc. Node/RAAN [rad]
                             temp(5) * pi/180, ... % Arg. Perigee [rad]
                             temp(6) * pi/180];    % True Anomaly [rad]
                         
        end

        
        % Get the True Anomaly [deg]. The True Anomaly is computed from the
        % Mean Anomaly [deg]         
        function [theta, E] = getTrueAnomaly(obj)
                        
            % Get the Mean Anomaly [deg] and convert in radians
            M0_ = obj.M0 * pi/180 ; % [rad]
            
            % Get the Eccentricity
            ecc = obj.e;
            
            % Compute the Eccentricity Anomaly
            E = CalcEA(M0_, ecc); % [rad]            
            
            % Compute the True Anomaly
            
            % NEW APPROACH (THAT IS, ORIGINAL ONE) - MARILENA 5th MAY
            cos_theta = ( cos(E) - ecc ) / (1 - ecc*cos(E));
            sin_theta = ( sin(E) * sqrt(1-ecc^2) ) / (1 - ecc*cos(E));
            
            theta = atan2(sin_theta, cos_theta);
                        
            % Convert the True Anomaly to Degrees
            theta = mod(theta *180/pi, 360); % in DEG

            % Convert the Eccentricity Anomaly in degrees
            E = E * 180/pi; % [deg]
             
        end % End GetTrueAnomaly() Method
            
        
    end % End Methods Section
    
    
end % End KeplerianElements Class

