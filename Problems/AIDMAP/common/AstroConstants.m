% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%-----------Copyright (C) 2018 University of Strathclyde and Authors-----------
%
%
%
%% AstroConstants: This class holds the astrodynamic constants
% 
% Inputs:
% *
% 
% Outputs: 
% * 
% 
%% Author(s): Juan Manuel Romero Martin (2014)
% Email:  juan.romero-martin@strath.ac.uk

classdef AstroConstants
    
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

