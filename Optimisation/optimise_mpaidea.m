function [x, fval, exitflag, output] = optimise_mpaidea(fitnessfcn, LB, UB, options)

%% optimise_mpaidea

%% Inputs:
%
% * fitnessfcn : function handle to function to optimise
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
[memories_record, memories, archivebest, population_evolution, vval_evolution,...
    B_mean, delta_local, options, exitflag] = MP_AIDEA_ALR(fitnessfcn, LB, UB, options.population, options);

% Output
x    = zeros(size(memories_record,3), numel(LB));
fval = zeros(size(memories_record,3), 1);

for i = 1 : size(memories_record,3)
    x(i,:)  = memories_record(end, 1:end-1, i);
    fval(i) = memories_record(end, end, i);
end

output.memories             = memories;
output.archivebest          = archivebest;
output.population_evolution = population_evolution;
output.vval_evolution       = vval_evolution;
output.B_mean               = B_mean;
output.delta_local          = delta_local;

end
