%------------------------------------------------------------------------%
% Copyright 2013 -  Juan Manuel Romero Martin
%------------------------------------------------------------------------%
%
% The Asteroids Class gather a collision of NEO ... blablabla 
%
% Author:
%   Juan Manuel Romero Martin
%
% Reference:
%   The Orbital elements of the Atira asteroids have been obtained from the
%   NEODys (Orb9 propagator) <http://newton.dm.unipi.it/neodys/>
%
% TODO:
%   - Finish the documentation of the code
%

classdef Asteroids 
    
     properties(Constant = true)  
         
        %Earth
        Root = CelestialBody('Root',    ... % Name 
                                  1, ... % Semimajor axis [AU]  
                                  0.0167, ... % Eccentricity 
                                  0,  ... % Inclination [deg]  
                                  -11.26064, ... % Asc. Node/raan [deg]
                                  102.94719, ... % Arg. Perigee [deg]
                                  124.9799993,   ... % Mean anomoly, M at time given t0 [deg]
                                  5972.5);          % Epoch [MJD2000]
         
        % Asteroid '163693'
        neo163693 = CelestialBody('163693',    ... % Name 
                                  0.741089104, ... % Semimajor axis [AU]  
                                  0.322141055, ... % Eccentricity 
                                  25.6173096,  ... % Inclination [deg]  
                                  103.9244961, ... % Asc. Node/raan [deg]
                                  252.9303171, ... % Arg. Perigee [deg]
                                  78.835491,   ... % Mean anomoly, M at time given t0 [deg]
                                  5055.5);          % Epoch [MJD2000]

        % Asteroid '164294'
        neo164294 = CelestialBody('164294',    ... % Name 
                                  0.617616364, ... % Semimajor axis [AU] 
                                  0.454500838, ... % Eccentricity 
                                  2.9495482,   ... % Inclination [deg]  
                                  211.4173031, ... % Asc. Node/raan [deg]
                                  5.1634356,   ... % Arg. Perigee [deg]
                                  290.9933747, ... % Mean anomoly, M at time given t0 [deg]
                                  5055.5);          % Epoch [MJD2000]

        % Asteroid '1998DK36' 
        neo1998DK36 = CelestialBody('1998DK36',  ... % Name  
                                    0.691925623, ... % Semimajor axis [AU] 
                                    0.415934335, ... % Eccentricity 
                                    2.025134,    ... % Inclination [deg]
                                    151.0774598, ... % Asc. Node/raan [deg]
                                    180.4707574, ... % Arg. Perigee [deg]
                                    271.3336628, ... % Mean anomoly, M at time given t0 [deg]
                                    5055.5);          % Epoch [MJD2000]

        % Asteroid '2004JG6' 
        neo2004JG6 = CelestialBody('2004JG6',   ... % Name 
                                   0.635211196, ... % Semimajor axis [kmn] 
                                   0.531143051, ... % Eccentricity 
                                   18.9451605,  ... % Inclination [deg]  
                                   37.0434558,  ... % Asc. Node/raan [deg]
                                   352.9835816, ... % Arg. Perigee [deg]
                                   109.9962091, ... % Mean anomoly, M at time given t0 [deg]
                                   5055.5);          % Epoch [MJD2000]

        % Asteroid '2005TG45' 
        neo2005TG45 = CelestialBody('2005TG45',  ... % Name 
                                    0.681458188, ... % Semimajor axis [AU] 
                                    0.372355361, ... % Eccentricity 
                                    23.3295341,  ... % Inclination [deg]  
                                    273.4637247, ... % Asc. Node/raan [deg]
                                    230.4110166, ... % Arg. Perigee [deg]
                                    344.1719563, ... % Mean anomoly, M at time given t0 [deg]
                                    5055.5);          % Epoch [MJD2000]

        % Asteroid '2006WE4'
        neo2006WE4 = CelestialBody('2006WE4',   ... % Name 
                                   0.784689884, ... % Semimajor axis [AU] 
                                   0.182947155, ... % Eccentricity 
                                   24.7676581,  ... % Inclination [deg]  
                                   311.038212,  ... % Asc. Node/raan [deg]
                                   318.6036423, ... % Arg. Perigee [deg]
                                   197.5386257, ... % Mean anomoly, M at time given t0 [deg]
                                   5055.5);          % Epoch [MJD2000]

        % Asteroid '2007EB26'
        neo2007EB26 = CelestialBody('2007EB26',  ... % Name 
                                    0.547815771, ... % Semimajor axis [AU] 
                                    0.788043512, ... % Eccentricity 
                                    8.4731244,   ... % Inclination [deg]  
                                    63.0561011,  ... % Asc. Node/raan [deg]
                                    236.8440102, ... % Arg. Perigee [deg]
                                    24.5936245,  ... % Mean anomoly, M at time given t0 [deg]
                                    5055.5);          % Epoch [MJD2000]

        % Asteroid '2008EA32'
        neo2008EA32= CelestialBody('2008EA32',  ... % Name 
                                   0.615950228, ... % Semimajor axis [AU] 
                                   0.304924579, ... % Eccentricity 
                                   28.2650298,  ... % Inclination [deg]  
                                   100.9714144, ... % Asc. Node/raan [deg]
                                   181.8480687, ... % Arg. Perigee [deg]
                                   120.9897376, ... % Mean anomoly, M at time given t0 [deg]
                                   5055.5);          % Epoch [MJD2000]

        % Asteroid '2008UL90'
        neo2008UL90 = CelestialBody('2008UL90',  ... % Name 
                                    0.694803229, ... % Semimajor axis [AU] 
                                    0.380248776, ... % Eccentricity 
                                    24.3094022,  ... % Inclination [deg]  
                                    81.1721556,  ... % Asc. Node/raan [deg]
                                    183.5991464, ... % Arg. Perigee [deg]
                                    37.7690364,  ... % Mean anomoly, M at time given t0 [deg]
                                    5055.5);          % Epoch [MJD2000]

        % Asteroid '2010XB11'
        neo2010XB11 = CelestialBody('2010XB11',  ... % Name 
                                    0.618033263, ... % Semimajor axis [AU] 
                                    0.533824731, ... % Eccentricity 
                                    29.8832556,  ... % Inclination [deg]  
                                    96.3248879,  ... % Asc. Node/raan [deg]
                                    202.473502,  ... % Arg. Perigee [deg]
                                    132.3162587, ... % Mean anomoly, M at time given t0 [deg]
                                    5055.5);          % Epoch [MJD2000]

        % Asteroid '2012VE46'
        neo2012VE46 = CelestialBody('2012VE46',  ... % Name 
                                    0.712824453, ... % Semimajor axis [AU] 
                                    0.361506763, ... % Eccentricity 
                                    6.6661913,   ... % Inclination [deg]  
                                    8.9484845,   ... % Asc. Node/raan [deg]
                                    190.3657028, ... % Arg. Perigee [deg]
                                    73.3979587,  ... % Mean anomoly, M at time given t0 [deg]
                                    5055.5);          % Epoch [MJD2000]

        % Asteroid '2013JX28'
        neo2013JX28 = CelestialBody('2013JX28',  ... % Name 
                                    0.600893413, ... % Semimajor axis [AU] 
                                    0.563999991, ... % Eccentricity 
                                    10.7666176,  ... % Inclination [deg]  
                                    39.9921724,  ... % Asc. Node/raan [deg]
                                    354.847281,  ... % Arg. Perigee [deg]
                                    185.3451515, ... % Mean anomoly, M at time given t0 [deg]
                                    5055.5);          % Epoch [MJD2000]
                                
                                
        %Additional Atira asteroids found                        
        neo2013TQ5 = CelestialBody('2013TQ5',  ... % Name 
                                    0.773678583, ... % Semimajor axis [AU] 
                                    0.155608989, ... % Eccentricity 
                                    16.39858808,  ... % Inclination [deg]  
                                    286.7788721,  ... % Asc. Node/raan [deg]
                                    247.3049075,  ... % Arg. Perigee [deg]
                                    232.5338344, ... % Mean anomoly, M at time given t0 [deg]
                                    6055.5);          % Epoch [MJD2000]
                                
        neo2014FO47 = CelestialBody('2014FO47',  ... % Name 
                                    0.75217103, ... % Semimajor axis [AU] 
                                    0.271103034, ... % Eccentricity 
                                    19.1979544,  ... % Inclination [deg]  
                                    358.659958,  ... % Asc. Node/raan [deg]
                                    347.4557791,  ... % Arg. Perigee [deg]
                                    52.10898445, ... % Mean anomoly, M at time given t0 [deg]
                                    6055.5);          % Epoch [MJD2000]
                                                                
        neo2015DR215 = CelestialBody('2015DR215',  ... % Name 
                                    0.666373934, ... % Semimajor axis [AU] 
                                    0.471605157, ... % Eccentricity 
                                    4.090341693,  ... % Inclination [deg]  
                                    314.9818619,  ... % Asc. Node/raan [deg]
                                    42.26044456,  ... % Arg. Perigee [deg]
                                    50.88874632, ... % Mean anomoly, M at time given t0 [deg]
                                    6055.5);          % Epoch [MJD2000]
                                                                                                
        neo2015ME131 = CelestialBody('2015ME131',  ... % Name 
                                    0.804887454, ... % Semimajor axis [AU] 
                                    0.198922545, ... % Eccentricity 
                                    28.87648991,  ... % Inclination [deg]  
                                    314.363762,  ... % Asc. Node/raan [deg]
                                    164.0284595,  ... % Arg. Perigee [deg]
                                    189.7431366, ... % Mean anomoly, M at time given t0 [deg]
                                    5652.5);          % Epoch [MJD2000]
                                
                                
                                
         
     end    
       
end

