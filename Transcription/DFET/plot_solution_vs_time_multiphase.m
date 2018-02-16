function plot_solution_vs_time_multiphase(x_in,problem,overlap,varargin)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%
% Plot entire solution of multiphase problem. Can either plot all states
% and control on same figure, or have a figure per phase.

if isempty(varargin)
       
    plotid = 1;
        
else

    if length(varargin)>1
   
        error('Only 1 extra parameter expected (this figure id)');
        
    else    
    
        if mod(varargin{1},1)>0
        
            error('Figure id must be an integer');        
    
        else
        
            plotid = varargin{1};
        
        end
    
    end
    
end

mint0 = inf;
maxtf = -inf;

minu = inf;
maxu = -inf;

figure

for i = 1:problem.num_phases
    
    %% Get normalised solution
    
    % Get chunk of solution relative to this phase
    
    x_sol = x_in(problem.phase_mask==i);
    x_sol = x_sol(:);                           % ensure column vector
    
    % Rescale problem
    x_sol = x_sol.*problem.structure{i}.scales.scale_opt;
    
    if problem.structure{i}.imposed_t0==1
        
        t_0 = problem.structure{i}.t_0;
        
    else
        
        t_0 = x_sol(problem.structure{i}.t0_vars);
        
    end
    
    if problem.structure{i}.imposed_tf==1
        
        t_f = problem.structure{i}.t_f;
        
    else
        
        t_f = x_sol(problem.structure{i}.tf_vars);
        
    end
    
    x0_best = problem.structure{i}.x_0;
    x0_best(~problem.structure{i}.imposed_initial_states) = x_sol(problem.structure{i}.x0_vars);
    els_best = problem.structure{i}.uniform_els;
    xx = x_sol(~problem.structure{i}.static_vars);
    [x_best,u_best,x_b] = extract_solution(xx,problem.structure{i},problem.structure{i}.x_f);
    
    %% Plot best solution
    
    if overlap == 1
    
        plot_solution_vs_time(x_best,u_best,x0_best,x_b,t_0,t_f,els_best,problem.structure{i},plotid);
        
    else
        
        plot_solution_vs_time(x_best,u_best,x0_best,x_b,t_0,t_f,els_best,problem.structure{i},plotid+i);

        
    end
    
    mint0 = min(mint0,t_0);
    maxtf = max(maxtf,t_f);
    
%     if ~isempty(problem.structure{i}.control_bounds)
% 
%             minu = min([minu min(problem.structure{i}.control_bounds(:,1))]);
%             maxu = max([maxu max(problem.structure{i}.control_bounds(:,2))]);
% 
%     end
    
    
    
end

% subplot(2,1,2)
% 
% if maxu<minu    % can happen with NO controls at all
% 
%     maxu= inf;
%     minu= -inf;
%     
% % else
% %     
% %     minu = minu*0.9;
% %     maxu = maxu/0.9;
% 
% end
% 
% axis([mint0 maxtf minu maxu]);

end