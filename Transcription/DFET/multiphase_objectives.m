function [val,dval] = multiphase_objectives (x_in,problem,jacflag)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%
% Computes objective function values and their gradient, for a multiphase
% problem.

% WE STILL DON'T HAVE FUNCTIONS OF GLOBAL STATIC VARIABLES NOT ASSOCIATED
% TO ANY PHASE. MAYBE LATER, BECAUSE WE COULD ALSO NEED TO PASS THAT GLOBAL
% STATIC VARIABLE DOWN TO THE PHASES, AND COULD BE A BIT MESSY

objectives = cell(problem.num_phases,1);        % create a cell array of all objectives of all phases
grad_objectives = cell(problem.num_phases,1);   % and another for their gradients

% Evaluate all objective function for each phase

for i = 1:problem.num_phases
    
    x_this_phase = x_in(problem.phase_mask==i); % get variables of this phase
    
    [objectives{i},grad_objectives{i}] = general_objectives (x_this_phase,problem.structure{i},jacflag);
    
end

% Call "outer level objective function" (function of the objective function
% of the individual phases, a generic function linking objectives of
% different phases)

[val,dval] = join_objectives(objectives,grad_objectives,jacflag,problem);

end

