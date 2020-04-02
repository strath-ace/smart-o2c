function [C,CEQ]=ps_constr(y,func, cfunc, lambda,z,zstar,scales,varargin)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%------------------------ Massimiliano Vasile -------------------------
%
% Pascoletti-Serafini scalarisation

%  INPUT
%  y=[t x]
%  func , cfunc - objective and constraint function handle
%  lambda - weights
%  z - reference point
%  zstar - normalisation point
%  arg - additional input to func
%
%  OUTPUT
%  C - constraint vector

% Massimiliano Vasile 2017, Lorenzo Ricciardi 2017

CEQ=0;
C=[];

t=y(1);
x=y(2:end).*scales;
f=func(x,varargin{:});
if ~isempty(cfunc)
    [C,CEQ]=cfunc(x,varargin{:});
end

f2 = (f-z)./(((zstar-z).*((zstar-z)~=0))+z.*(abs(z)>=1).*((zstar-z)==0)+(abs(z)<1).*((zstar-z)==0));

tv = ones(size(lambda))*t;
tv(lambda==-1) = 0;
lambda = abs(lambda);

C=[tv(1)-max(lambda);(lambda.*f2-tv)';C];