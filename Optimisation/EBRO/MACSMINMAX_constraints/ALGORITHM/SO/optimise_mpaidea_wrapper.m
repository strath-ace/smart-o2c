% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [ x, fval, exitflag, output ] = optimise_mpaidea_wrapper(problem,par)

if (isfield(par,'initial_population') && ~isempty(par.initial_population))
    initial_population = par.initial_population;
else
    % Initialise populations
    initial_population = zeros(par.n_agents,problem.dim,par.n_populations);

    for s = 1 : par.n_populations
        % pop = lhsdesign(par.n_agents,problem.dim,'criterion','maximin').*repmat(problem.ub-problem.lb,par.n_agents,1)+repmat(problem.lb,par.n_agents,1);
        pop = lhsgen(par.n_agents,problem.dim).*repmat(problem.ub-problem.lb,par.n_agents,1)+repmat(problem.lb,par.n_agents,1);

        
        initial_population(:,:,s) = pop;
    end
end




fitnessfcn = problem.fitnessfcn;



% Run MP-AIDEA
[memories_record, memories, archivebest, archiveALL, population_evolution, vval_evolution,...
    B_mean, delta_local, inite, iglob, options, exitflag] = MP_AIDEA(fitnessfcn, problem.lb, problem.ub, initial_population, par, problem.par_objfun);



% Output: minima and minima's objective value
x    = zeros(size(memories_record,3), problem.dim);
fval = zeros(size(memories_record,3), 1);

for i = 1 : size(memories_record,3)
    x(i,:)  = memories_record(end, 1:end-1, i);
    fval(i) = memories_record(end, end, i);
end



output.memories_record      = memories_record;
output.memories             = memories;
output.archivebest          = archivebest;
output.archiveALL           = archiveALL;
output.population_evolution = population_evolution;
output.vval_evolution       = vval_evolution;
output.B_mean               = B_mean;
output.delta_local          = delta_local;
output.number_LR            = inite;
output.number_GR            = iglob;
output.nfeval               = par.nFeValMax;
return