% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%-----------Copyright (C) 2016 University of Strathclyde-------------
%
%
%
% test_cec09.m
% 
% Matlab source codes
% 
% Test function evaluation for CEC 2009 MOO Competition
%
% Usage: test_cec09(problem_name, vairaible_dimension), for example
% test_cec09('CF1',10)
% 
% Please refer to the report for more information.
%

function [f]=test_func_cec09(xi,name, dim)

% get the search boundary
xrange = xboundary(name, dim);

% randmoly generate NN points in the search space
x  = xrange(:,1) + (xrange(:,2)-xrange(:,1)).*xi';

% %% === test matlab version with constraints
% % get the function handle
% fobj = cec09(name);
% % test the function
% [y,c] = fobj(x);
% % display the results
% disp('objectives');  disp(y);
% disp('constraints'); disp(c);

%% === test matlab version withot constraints
% get the function handle
fobj = cec09(name);
% test the function
[f] = fobj(x);
f=f';
end
