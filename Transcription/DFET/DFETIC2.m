function x_guess = DFETIC2(x_0,t_0,t_f,static,u,structure)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%
% Construction of initial guess by DFET, all elements simultaneously.

% The simplest way to do it is RE-transcribing part of the
% problem, might be inefficient but it's cleaner and more generic...

%% Clone original structure

temp_structure = structure;
temp_structure.constants = structure.constants;

%% generate initial guess (cloning IC)

temp_structure.init_type = 'CloneIC';
x_guess = make_first_guess(x_0,t_0,t_f,static,u,temp_structure);

x_guess = x_guess./temp_structure.scales.scale_opt;

% Check variables are within bounds. They should never be outside at this
% stage

if any(x_guess<=temp_structure.norm_lbv) || any(x_guess>=temp_structure.norm_ubv)
   
    keyboard
    
end

% Solve this feasibility problem

options = optimoptions(@fmincon,'Algorithm','interior-point','Display','iter-detailed','MaxFunEvals',10*length(x_guess),'MaxIter',10*length(x_guess),'GradConstr','on','GradObj','on','TolFun',1e-6,'Tolcon',1e-6,'TolX',1e-6,'DerivativeCheck','off');

doit = 1;
maxit = 1;
it = 0;

while doit && (it<maxit)

    x_guess = fmincon(@(x) feas_only(x),x_guess,[],[],[],[],temp_structure.norm_lbv,temp_structure.norm_ubv,@(x) Full_constr_dynamics(temp_structure,x,1),options);                                                                                                                                                                             
        
    [c,ceq] = Full_constr_dynamics(temp_structure,x_guess,0);
    newfeas = max([abs(ceq);c]);
    it = it+1;
    
    if (newfeas<=options.TolCon)
        
        doit = 0;
       
    end
        
end

x_guess = x_guess.*structure.scales.scale_opt;

end
