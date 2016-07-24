%% StardustTool: This class contains various functions related  to finding the Keplerian and Cartesian elements at a certain time for a certain orbit
% 
%% Inputs:
% * 
% 
%% Outputs: 
% * 
% 
%% Author(s): Juan Manuel Romero Martin (2014), Marilena Di Carlo (2014)
% Email:  juan.romero-martin@strath.ac.uk, marilena.di-carlo@strath.ac.uk

classdef StardustTool
    
    methods (Static = true)
        
        % Get the Orbital Elements at certain specific EPOCH (MJD2000)
        % 
        % INPUT:
        %  celestial_obj  The Celestial Body   
        %  epoch          The Epoch [MDJ2000]
        % 
        % OUTPUT:
        % 
        % Keplerian Elements Object:
        %    a   Semimajor axis [AU]     
        %    e   Eccentricity 
        %    i   Inclination [deg]     
        %    OM  Asc. Node/raan [deg]
        %    W   Arg. Perigee [deg]     
        %    M   Mean anomoly, M at time given Epoch [deg]
        %    t   Time at which M is computed [MJD] 

        function kep = KeplerianElementsAt2(celestial_obj, epochMJD2000)
         
            % Sun planetary constant [km^3/s^2]
            mu_sun = AstroConstants.Sun_Planetary_Const;            

            % Point to the Orbital Elements of the object (CLEAN THIS)
            a_  = celestial_obj.a *  AstroConstants.Astronomical_Unit; % [km]
            e_  = celestial_obj.e;
            di_ = celestial_obj.i;
            omm_= celestial_obj.OM;
            om_ = celestial_obj.W;
            m_  = celestial_obj.M0;
            
            % compute the Mean Motion [1/sec]
            n_  = sqrt(mu_sun/a_^3);
                                   
            % Compute the different time [MDJ2000] --> [sec]
            diffTime = (epochMJD2000 - celestial_obj.t0) * 86400;
            
            % Mean Anomaly at the given time [rad]
            Mf = m_ * pi/180 + n_ * diffTime;
            
            Mf = mod(Mf * 180/pi, 360);
            
            % Save the Semi-major in AU
            a_  = celestial_obj.a;  % [AU]
                        
            % Create the Output Keplerian Element
            kep = KeplerianElements(a_,   ...      % Semimajor axis [AU] 
                                    e_,   ...      % Eccentricity 
                                    di_,  ...      % Inclination [deg]  
                                    omm_, ...      % Asc. Node/raan [deg]
                                    om_,  ...      % Arg. Perigee [deg]     
                                    Mf, ...  % Mean anomoly, M at time given epoch [deg]                                                               
                                    epochMJD2000); % Epoch [MJD2000]
            
        end
                                
        
        % Get the Cart
        % 
        % INPUTS
        % 
        % OUTPUTS
        % R Position Vector [km]
        % V Velocity Vector [km/s]
        % 
        % TODO
        % - Change the name of the method. 
        % 
        function [r, v] = CartesianElementsAt(celestial_obj, epochMDJ2000)
                            
            % Compute the Cartesian Vectors for Celestial Body Object
            if isa(celestial_obj, 'CelestialBody')

                % Compute the Keplerian Elements at given Epoch
                kep = StardustTool.KeplerianElementsAt2(celestial_obj, epochMDJ2000);
                                                  
                % Transform from Kepler Orbital Elements to cartesian                                              
                akep = kep.getKepArray_KM_RAD();

                % Convert From Keplerian to Cartesian
                car = kep2cart(akep, ...                            % Kepler Orbital Elements [km and rad]
                               AstroConstants.Sun_Planetary_Const); % Gravitational constant of the Sun [km^3/s^2]
                       
                r   = car(1:3); % [km]
                v   = car(4:6); % [km/s]
                
            % Compute the Cartesian Vectors for the Planets                
            elseif isa(celestial_obj, 'double')
                
                % Get teh Peh
                [r, v] = EphSS( celestial_obj, epochMDJ2000 );
                
            end
                
        end % End CartesianElementsAt Method     
               
            
    end % End STATIC Methods
    
end % End of StardustTool Class

