function [c,ceq,Jc,Jceq] = Pascoletti_Serafini_constraints(in,lambda,z,zstar,options,tstart)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%
% Pascoletti-Serafini scalarisation

alpha_0 = max(lambda);
alpha = in(1);
x = in(2:end);

% Eval constraints and Jacobians wrt dynamic variables and time (NORMALISED
% VARIABLES)

[c,ceq,Jc,Jceq] = multiphase_constr_full(options.oc.problem,x,1,tstart,options.timeout_single_level_opt);

Jceq = [zeros(1,length(ceq)); Jceq];  % zeros is the derivative of equality constraints wrt alpha

% Eval objectives 

[cc,Jcc] = multiphase_objectives(x,options.oc.problem,1);

cc = cc';   % cc is returned as a column vector, but we need a row vector here

% cast objectives as constraints (i.e. normalised Pascoletti-Serafini scalarisation)

alphav = ones(size(lambda))*alpha;

alphav(lambda==-1) = 0;


vvv = -ones(size(lambda));
vvv(lambda==-1) = 0;

lambda = abs(lambda);

cPS = (lambda.*(cc-z)./(((zstar-z).*((zstar-z)~=0))+((zstar-z)==0))-alphav)';

qq = repmat(lambda./(((zstar-z).*((zstar-z)~=0))+((zstar-z)==0)),size(Jcc,1),1); % multiplier lambda/(z*-z) to get correct Jacobian of constraints
Jc = [vvv zeros(1,size(Jc,2)); qq.*Jcc Jc]; % add derivatives wrt alpha and scale for normalised objectives

c = [alpha-alpha_0; cPS; c]; % alpha must be lower than initial value for a feasible solution

v = zeros(size(Jc,1),1); % derivative of alpha-alpha_0 wrt all variables (i.e. [1 0 0 0 0... 0]
 
v(1) = 1;
 
Jc = [v Jc];        % add gradient of alpha-alpha_0 constraint

end