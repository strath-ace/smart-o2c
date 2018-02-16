function [xt,ut,t] = eval_solution_over_time_multiphase(x_in,problem,denormalise)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%
% Evaluates the solution on a time interval for multiphase problems

[x,u] = extract_solution_multiphase(x_in,problem,denormalise);

xt = cell(problem.num_phases,1);
ut = cell(problem.num_phases,1);
t =  cell(problem.num_phases,1);

for i = 1:problem.num_phases
    
    xx = x_in(problem.phase_mask==i);  
    
    if denormalise 
       
        xx = xx.*problem.structure{i}.scales.scale_opt;
        
    end
    
    % get t0
    
    if problem.structure{i}.imposed_t0                
        
        t_0 = problem.structure{i}.t_0;
        
    else
        
        t_0 = xx(problem.structure{i}.t0_vars);            
        
    end               
        
    % get t_f
    
    if problem.structure{i}.imposed_tf
        
        t_f = problem.structure{i}.t_f;
        
    else
        
        t_f = xx(problem.structure{i}.tf_vars);
            
    end
    
    % create t
    qq = problem.structure{i}.uniform_in_nodes_state*(t_f-t_0)+t_0;
    t{i} = sort(qq(:));
    %t{i} = linspace(t_0,t_f,1000);
    [xt{i},ut{i}] = eval_solution_over_time(x{i},u{i},t_0,t_f,t{i},problem.structure{i}.uniform_els,problem.structure{i});

end

end