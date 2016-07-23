function equ = kep2eq( kep )

% Convertion from keplerian elements to equinoctial elements.
% All the units to be consistent, angles in radians.
%
% INPUT
%   kep     Vector of Keplerian elements:
%               kep = [a, e, incl, Om, om, th]
%
% OUTPUT
%  equ     Vector of equinoctial elements:


a    = kep(1);
e    = kep(2);
incl = kep(3);
Om   = kep(4);
om   = kep(5);
th   = kep(6);


P1 = e * sin(Om+om);
P2 = e * cos(Om+om);
Q1 = tan(incl/2) * sin(Om);
Q2 = tan(incl/2) * cos(Om);
L = mod(Om + om + th, 2*pi);


equ = [a P1 P2 Q1 Q2 L];

end

