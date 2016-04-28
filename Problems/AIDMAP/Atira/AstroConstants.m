classdef AstroConstants
    
    % ASTROCONSTANTS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Constant = true)  
        
        % Conversion Factor from DEG to RAD
        DEG2RAD = pi / 180;
        
        % Conversion Factor from RAD to DEG
        RAD2DEG = 180 / pi;
        
        % Universal gravity constant (G) [km^3/(kg*s^2)]
        Universal_Gravity_Const = 6.67259e-20; % From DITAN
        
        % Astronomical Unit (AU) [km]
        Astronomical_Unit = 149597870.7; % From DITAN
        
        % Sun mean radius [km]
        Sun_Mean_Radius = 700000; % From DITAN
      
        % Sun Gravitational parameter (mu = mass * G) [km^3/s^2]
        Sun_Planetary_Const = 0.19891000000000E+31*6.67259e-20; % From DITAN
        
    end
    
    methods        
        % NO METHODS ARE REQUIRED    
    end
    
end

