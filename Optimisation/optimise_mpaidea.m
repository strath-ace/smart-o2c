function [x,fval,exitflag,output] = optimise_mpaidea(fitnessfcn, nvars, LB, UB, options)

%% optimise_mpaidea

%% Inputs:
%
% * fitnessfcn : function handle to function to optimise
% * nvars : number of dimension of the problem
% * LB : lower boundaries
% * UB : upper boundaries
% * options : 
%
%% Output:
% * y : explanation
%
% Author: Marilena Di Carlo
% email: marilena.di-carlo@strath.ac.uk

% Check dimension
if length(LB) ~= length(UB) || length(LB) ~= nvars || length(UB) ~= nvars
    error('Dimension of upper or lower boundary not compatible with dimension of the problem')
end


% Run MP-AIDEA
[memories, B_mean, bubble, archivebest,options, exitflag] = MP_AIDEA_ALR(fitnessfcn, LB, UB, options.population, options);

% Define output
x = memories(1,1:nvars);
fval = memories(1,nvars + 1);
output = [];
end
