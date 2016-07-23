function h = plotorbit(kep5, mu, s)

% Plot of the orbit
%
%   h = plotorbit(kep5, mu, s)
%
% It plots a keplerian orbit without integrating the trajectory, but using
% [a, e, i, Om, om]
%
% INPUT:
%       kep5 = 5 keplerian parameters of the orbit [a, e, i, Om, om]
%              (angles in rad)
%       mu = planetary constant (mu = mass * G) [L^3/T^2]
%       s = character string made from one element from any or all the 
%           following 3 columns:
%               b     blue          .     point              -     solid
%               g     green         o     circle             :     dotted
%               r     red           x     x-mark             -.    dashdot 
%               c     cyan          +     plus               --    dashed   
%               m     magenta       *     star             (none)  no line
%               y     yellow        s     square
%               k     black         d     diamond
%                                   v     triangle (down)
%                                   ^     triangle (up)
%                                   <     triangle (left)
%                                   >     triangle (right)
%                                   p     pentagram
%                                   h     hexagram
%
% OUTPUT:
%       h = handle to the plot object.
%
% functions called: kep2cart
%
% - Camilla Colombo - 17/02/2006
%                   - 12/02/2007 - changed kep2cart prototype
%
% - Revised by Matteo Ceriotti - 20/11/2007
% - Matteo Ceriotti - 20/11/2007 - Added handle to plot object as output
%                                - Added initialisation of cart matrix for
%                                  speed
%
% ------------------------- - SpaceART Toolbox - --------------------------

if nargin < 3
    s = 'b';
end

kep5 = kep5(:);
kep5 = kep5(1:5);

rdeg = pi/180;
cart = zeros(360, 6);
for i = 1:360
    f = i*rdeg;
    cart(i, :) = kep2cart([kep5;f], mu);
end
h = plot3(cart(:, 1), cart(:, 2), cart(:, 3), s);

xlabel('x');
ylabel('y');
zlabel('z');

return