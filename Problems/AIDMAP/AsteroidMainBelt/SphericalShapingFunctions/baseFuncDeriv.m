function [R,Phi,RPrime1,PhiPrime1,RPrime2,PhiPrime2,RPrime3,PhiPrime3,D,DPrime,isDpos] = baseFuncDeriv(theta,param)
% Function that calculates the functions R and Phi and their derivatives
% until the third with the angle theta and the parameters
% INPUT: angle theta (can be vector), vector param containing the parameters in the
% following order: param = x + a2 = [a0 a1 a3..a6 b0..b3 a2]

% OUTPUT: R,Phi and their derivatives until the third, functions D and D'
% isDpos: 1 if D>0 for every theta, 0 otherwise. They are vectors if theta is
% a vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Niccolo' Gastaldello, October 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute
isDpos = 1;

% Parameters assignment
a0 = param(1);
a1 = param(2);
a2 = param(end);
a3 = param(3);
a4 = param(4);
a5 = param(5);
a6 = param(6);
b0 = param(7);
b1 = param(8);
b2 = param(9);
b3 = param(10);

% Compute R,Phi and their derivatives
R = 1./(a0 + a1*theta + a2*theta.^2 + (a3 + a4*theta).*cos(theta) + (a5 + a6*theta).*sin(theta));
Phi = (b0+b1*theta).*cos(theta) + (b2+b3*theta).*sin(theta);

RPrime1 = -R.^2.*(a1 + 2*a2*theta - a3*sin(theta) + a4*(cos(theta) - theta.*sin(theta)) + a5*cos(theta) + a6*(sin(theta) + theta.*cos(theta)));
PhiPrime1 = -b0*sin(theta) + b1*(cos(theta) - theta.*sin(theta)) + b2*cos(theta) + b3*(sin(theta) + theta.*cos(theta));

RPrime2 = -R.^2.*(2*a2 - a3*cos(theta) + a4*(-2*sin(theta) - theta.*cos(theta)) - a5*sin(theta) + a6*(2*cos(theta) - theta.*sin(theta))) + 2*(RPrime1.^2)./R;
PhiPrime2 = -b0*cos(theta) + b1*(-2*sin(theta) - theta.*cos(theta)) - b2*sin(theta) + b3*(2*cos(theta) - theta.*sin(theta));

RPrime3 = 6*RPrime1.^3./R.^2 - 6*R.*RPrime1.*(2*a2 - a3*cos(theta) + a4*(-2*sin(theta) - theta.*cos(theta)) - a5*sin(theta) + a6*(2*cos(theta) - theta.*sin(theta))) +...
             - R.^2.*(a3*sin(theta) + a4*(-3*cos(theta) + theta.*sin(theta)) - a5*cos(theta) + a6*(-3*sin(theta) - theta.*cos(theta)));
PhiPrime3 = b0*sin(theta) + b1*(-3*cos(theta) + theta.*sin(theta)) - b2*cos(theta) + b3*(-3*sin(theta) - theta.*cos(theta));

fact = PhiPrime1.^2 + cos(Phi).^2;
D = -RPrime2 + 2*(RPrime1.^2)./R + RPrime1.*PhiPrime1.*(PhiPrime2 - sin(Phi).*cos(Phi))./fact + R.*fact;
DPrime = -RPrime3 + 4*(RPrime1.*RPrime2)./R - 2*(RPrime1.^3)./R.^2 + ...
          +(RPrime2.*PhiPrime1.*(PhiPrime2 - sin(Phi).*cos(Phi)) + RPrime1.*PhiPrime2.*(PhiPrime2 - sin(Phi).*cos(Phi)) + RPrime1.*PhiPrime1.*(PhiPrime3 - PhiPrime1.*cos(Phi).^2 + PhiPrime1.*sin(Phi).^2) )./fact + ...
          -(RPrime1.*PhiPrime1.*(PhiPrime2 - sin(Phi).*cos(Phi)).*(2*PhiPrime1.*PhiPrime2 - 2*cos(Phi).*sin(Phi).*PhiPrime1))./fact.^2 + ...
          +RPrime1.*(fact) + R.*(2*PhiPrime1.*PhiPrime2 - 2*cos(Phi).*sin(Phi).*PhiPrime1);

% Check the curvature of the trajectory (with D) and discard the
% trajectory if D<0 for at least one angle (D depends only on the geometry
% of the trajectory)
for i=1:length(theta)
   if D(i) < 0
      isDpos = 0;
      break;
   end
end

end