function [val,dval] = general_objectives (x_in,structure,jacflag)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%
% Computes objective function value and its gradient

els = structure.uniform_els;
static = x_in(structure.other_vars);
nstatic = length(static);

if structure.imposed_t0
    
    t_0 = structure.t_0_norm;
    
else
    
    t_0 = x_in(structure.t0_vars);
    
end

if structure.imposed_tf
    
    t_f = structure.t_f_norm;
    
else
    
    t_f = x_in(structure.tf_vars);
    
end

n_x0 = sum(structure.x0_vars);  % number of free initial conditions
ids = 1:length(structure.x_0);
ids = ids(~structure.imposed_initial_states);   % get which x_0 are coded in solution vector
x_0 = structure.x_0_norm;
x_0(ids) = x_in(structure.x0_vars);

x_f = structure.x_f_norm;

x_sol = x_in(~structure.static_vars);

h = 1e-7;

if jacflag ==0
    
    dval = [];
    val = eval_cost_functions2(structure.g,structure.weights,x_sol,x_0,x_f,[t_0 t_f],static,els,structure,0,[],[],[],[]);  

else
    
    % value + derivatives wrt states, controls and final conditions
        
    [val,grad] = eval_cost_functions2(structure.g,structure.weights,x_sol,x_0,x_f,[t_0 t_f],static,els,structure,1,structure.dgu0,structure.dgxf,structure.dguf,structure.dgxi,structure.dgui);
    dval = zeros(length(x_in),length(val));
    dval(~structure.static_vars,:) = grad;  
    
    % derivatives wrt generic static variables
    
    if nstatic>0
        
        Jstatic = zeros(nstatic,length(val));
        
        for i=1:nstatic
            
            temp_static = static;
            temp_static(i) = temp_static(i)+h;
            f2 = eval_cost_functions2(structure.g,structure.weights,x_sol,x_0,x_f,[t_0 t_f],temp_static,els,structure,0,[],[],[],[]);
            Jstatic(i,:) = (f2-val)'/h;
            
        end

        dval(structure.other_vars,:) = Jstatic;
        
    end
   
    % derivative wrt initial time
    
    if any(structure.t0_vars)
        
        f2 = eval_cost_functions2(structure.g,structure.weights,x_sol,x_0,x_f,[t_0+h t_f],static,els,structure,0,[],[],[],[]);
        dfdt0 = (f2-val)'/h;
        dval(structure.t0_vars,:) = dfdt0;
        
    end
    
    % derivative wrt final time
    
    if any(structure.tf_vars)
        
        f2 = eval_cost_functions2(structure.g,structure.weights,x_sol,x_0,x_f,[t_0 t_f+h],static,els,structure,0,[],[],[],[]);
        dfdtf = (f2-val)'/h;
        dval(structure.tf_vars,:) = dfdtf;
        
    end
    
    % derivative wrt free x_0 
            
    if n_x0>0
        
        dfdx0 = zeros(n_x0,length(val));
        
        for i=1:n_x0
            
            temp_x0 = x_0;
            temp_x0(ids(i)) = temp_x0(ids(i))+h;
            f2 = eval_cost_functions2(structure.g,structure.weights,x_sol,temp_x0,x_f,[t_0 t_f],static,els,structure,0,[],[],[],[]);
            dfdx0(i,:) = (f2-val)'/h;
            
        end
        
        dval(structure.x0_vars,:) = dfdx0;

    end
            
end

end