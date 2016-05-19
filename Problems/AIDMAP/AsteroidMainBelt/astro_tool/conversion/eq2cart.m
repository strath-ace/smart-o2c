function cart = eq2cart (equin,j)

global mu

a  = equin(1);
P1 = equin(2);
P2 = equin(3);
Q1 = equin(4);
Q2 = equin(5);
L_E  = equin(6);
% 
% % Eccentricity
% e     = sqrt(P1^2 + P2^2);
% 
% % Parameter p
% p = a * (1-e^2);
% 
% alpha2 = Q1^2 - Q2^2;
% s2    = 1 + Q1^2 + Q2^2;
% w     = 1 + P1 * cos(L) + P2 * sin(L);
% r     = p/w;
% 
% cart(1) = (r/s2) * (cos(L) + alpha2 * cos(L) + 2 * Q1 * Q2 * sin(L));
% cart(2) = (r/s2) * (sin(L) - alpha2 * sin(L) + 2 * Q1 * Q2 * cos(L));
% cart(3) = (2*r/s2) * (Q1 * sin(L) - Q2 * cos(L));
% 
% cart(4) = -(1/s2) * (sqrt(mu/p)) * (sin(L)  + alpha2 * sin(L) - 2 * Q1 * Q2 * cos(L) + P2 - 2 * P1 * Q1 * Q2 + alpha2 * P2); 
% cart(5) = -(1/s2) * (sqrt(mu/p)) * (-cos(L) + alpha2 * cos(L) + 2 * Q1 * Q2 * sin(L) - P1 + 2 * P1 * Q1 * Q2 + alpha2 * P1); 
% cart(6) =  (2/s2) * (sqrt(mu/p)) * (Q1 * cos(L) + Q2 * sin(L) + P1 * Q1 + P2 * Q2); 

beta = 1 / (1 + sqrt(1 - P1^2 - P2^2) );

X1 = a * ( (1 - beta * P1^2) * cos(L_E) + P1 * P2 * beta * sin(L_E) - P2 );
Y1 = a * ( (1 - beta * P2^2) * sin(L_E) + P1 * P2 * beta * cos(L_E) - P1 );

n = sqrt( mu / a^3);
r = a * (1 - P2 * cos(L_E) - P1 * sin(L_E));

X1_dot = (n * a^2 / r) * ( beta * P1 * P2 * cos(L_E) - (1 - beta * P1^2) * sin(L_E) );
Y1_dot = (n * a^2 / r) * ( - beta * P1 * P2 * sin(L_E) + (1 - beta * P2^2) * cos(L_E) );

f = ( 1 / (1 + Q1 ^ 2 + Q2^2) ) * [1 - Q1^2 + Q2^2,  2 * Q1 * Q2,   -2 * Q1 * j];
g = ( 1 / (1 + Q1 ^ 2 + Q2^2) ) * [2 * Q1 * Q2 * j,      (1 + Q1^2 - Q2^2)*j,   2 * Q2];

cart(1:3) = X1 * f + Y1 *g;
cart(4:6) = X1_dot * f + Y1_dot *g;

end