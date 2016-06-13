function [TPrime1,TPrime2] = timeFunc(R,D,RPrime,DPrime)
% Function that calculates the time derivatives (first and second) for all
% the nodes of the trajectory computed at the precedent step
% INPUT: R,D and their first derivative
% OUTPUT: T' and T''

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Niccolo' Gastaldello, October 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Computation

global mu_S

TPrime1 = sqrt(D.*R.^2/mu_S);
TPrime2 = (2*D.*RPrime + R.*DPrime)./(2*sqrt(mu_S*D));
end

