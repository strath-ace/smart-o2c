function [x,fval,exitflag,output] = optimise_mpaidea(fitnessfcn, LB, UB, options)

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



% Run MP-AIDEA
[memories_out, memories, vval, B_mean, bubble, archivebest,...
    population_evolution, vval_evolution, options] = MP_AIDEA_ALR(fitnessfcn, LB, UB, options.population, options);

% Define output
% x = memories(1,1:nvars);
% fval = memories(1,nvars + 1);
% output = [];

keyboard
end
