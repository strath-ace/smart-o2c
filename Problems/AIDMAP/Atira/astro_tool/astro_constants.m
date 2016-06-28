function out=astro_constants(in)
% astro_constants - Returns astrodynamic-related physical constants
% 
%   out = astro_constants(in)
%
% Returns a row vector of constants, in which there is the corresponding
% constant for each element of the input vector.
%
% List of identifiers:
%   Generic astronomical constants:
%       1   Universal gravity constant (G) [km^3/(kg*s^2)]
%       2   Astronomical Unit (AU) [km]
%   Sun related:
%       3   Sun mean radius [km]
%       4   Sun planetary constant (mu = mass * G) [km^3/s^2]
%      31  Energy flux density of the Sun [W/m2 at 1 AU]
%   Other:
%       5   Speed of light in the vacuum [km/s]
%       6   Standard free fall (the acceleration due to gravity on the
%           Earth's surface at sea level) [m/s^2]
%       7   Mean distance Earth-Moon [km]
%       8   Obliquity (angle) of the ecliptic at Epoch 2000 [rad]
%   Planetary constants of the planets (mu = mass * G) [km^3/s^2]:
%       11  Me
%       12  V
%       13  E
%       14  Ma
%       15  J
%       16  S
%       17  U
%       18  N
%       19  P
%       20  Moon
%   Mean radius of the planets [km]:
%       21  Me
%       22  V
%       23  E
%       24  Ma
%       25  J
%       26  S
%       27  U
%       28  N
%       29  P
%       30  Moon
%
% INPUT
%   in      Vector of identifiers of required constants.
%
% OUTPUT
%   out     Vector of constants.
%
% EXAMPLES
%
%   astro_constants([2, 4, 26])
%      Returns a row vector in which there is the value of the AU, the Sun
%      planetary constant and the mean radius of Saturn.
%
%   astro_constants(10 + [1:9])
%      Returns a row vector with the planetary constant of each planet.
%
% Matteo Ceriotti, 2006                     Ver. 1.2
% Camilla Colombo, 26/10/06: updated
%                  22/10/07: astro_constants(8) added


% Notes for upgrading this function:
% 1) It is possible to add new constants. Please DO NOT change the
%    structure of this function, as well as its prototype.
% 2) DO NOT change the identifiers of the constants that have already been
%    defined in this function. If you want to add a new constant, please
%    use an unused identifier.
% 3) DO NOT add constants that can be easily computed starting form other
%    ones (avoid redundancy).
% 4) ALWAYS check that you are modifying the latest version of the
%    function, and make the new version available to everyone as soon as
%    possible.
% 5) When adding new constants, please upgrade the help of the function
%    and its version accordingly.
% Thanks!
%                                                       Matteo Ceriotti

for i=1:length(in)
    switch in(i)
        case 1
            out(i)=6.67259e-20; % From DITAN
        case 2
            out(i)=149597870.7; % From DITAN
        case 3
            out(i)=700000; % From DITAN
        case 4
            out(i)=0.19891000000000E+31*6.67259e-20; % From DITAN
        case 5
            out(i)=299792.458; % Definition in the SI
        case 6
            out(i)=9.80665; % Definition in Wertz
        case 7
            out(i)=384401; % Definition in Wertz
        case 8
            out(i)=23.43928111*pi/180; % Definition in Wertz
        case 11
            out(i)=0.33020000000000E+24*6.67259e-20; % From DITAN
        case 12
            out(i)=0.48685000000000E+25*6.67259e-20; % From DITAN
        case 13
            out(i)=0.59736990612667E+25*6.67259e-20; % From DITAN
        case 14
            out(i)=0.64184999247389E+24*6.67259e-20; % From DITAN
        case 15
            out(i)=0.18986000000000E+28*6.67259e-20; % From DITAN
        case 16
            out(i)=0.56846000000000E+27*6.67259e-20; % From DITAN
        case 17
            out(i)=0.86832000000000E+26*6.67259e-20; % From DITAN
        case 18
            out(i)=0.10243000000000E+27*6.67259e-20; % From DITAN
        case 19
            out(i)=0.14120000000000E+23*6.67259e-20; % From DITAN
        case 20
            out(i)=0.73476418263373E+23*6.67259e-20; % From DITAN
        case 21
            out(i)=0.24400000000000E+04; % From DITAN
        case 22
            out(i)=0.60518000000000E+04; % From DITAN
        case 23
            out(i)=0.63781600000000E+04; % From DITAN
        case 24
            out(i)=0.33899200000000E+04; % From DITAN
        case 25
            out(i)=0.69911000000000E+05; % From DITAN
        case 26
            out(i)=0.58232000000000E+05; % From DITAN
        case 27
            out(i)=0.25362000000000E+05; % From DITAN
        case 28
            out(i)=0.24624000000000E+05; % From DITAN
        case 29
            out(i)=0.11510000000000E+04; % From DITAN
        case 30
            out(i)=0.17380000000000E+04; % From DITAN
        case 31
            out(i)=1367; % From SMAD/Wertz
        % Add an identifier and constant here. Prototype:
        % case $identifier$
        %     out(i)=$constant_value$;
        otherwise
            warning('Constant identifier %d is not defined!',in(i));
            out(i)=0;
    end
end
