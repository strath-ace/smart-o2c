function [x,u,xb] = extract_solution_multiphase(x_in,problem,denormalise)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%
% From the vector solution of a multiphase problem, generates a cell vector
% for each phase. Each cell vector contains the nodal solution vector x,
% the nodal control vector u, and the boundary states vector xb

x = cell(problem.num_phases,1);
u = cell(problem.num_phases,1);
xb = cell(problem.num_phases,1);

for i = 1:problem.num_phases
    
    x_sol = x_in(problem.phase_mask==i);
    
    if denormalise
        
        x_sol = x_sol.*problem.structure{i}.scales.scale_opt; 
        
    end     
    
    x_sol = x_sol(~problem.structure{i}.static_vars);
    [x{i},u{i},xb{i}] = extract_solution(x_sol,problem.structure{i},problem.structure{i}.x_f);

end