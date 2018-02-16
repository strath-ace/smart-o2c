function [c,ceq,Jc,Jceq] = Full_constr_dynamics (structure,x_in,jacflag)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%
% This function computes all equality and inequality constraints: algebraic
% and differential path constraints, boundary points constraints, and their
% jacobians with respect to all states, controls and static variables.

[x_0,x_f,t_0,t_f,static,ids] = get_boundary_and_static_vars(x_in,structure);

n_x0 = sum(structure.x0_vars);  % number of free initial conditions
nstatic = length(static);       % number of static parameters
    
x_in = x_in(~structure.static_vars);
    
h = 1e-7;

%% Inequality constraints

c  = [];
Jc = [];

% Algebraic path constraints

if ~isempty(structure.c)
    
    [c,Jc] = eval_path_constraints(structure,x_in,x_0,x_f,t_0,t_f,static,structure.uniform_els,jacflag);
    
    if jacflag
        
        Jc = Jc';
        
        % derivatives wrt generic static variables
        
        if nstatic>0
            
            Jstatic = zeros(nstatic,size(Jc,2));
            
            for i=1:nstatic
                
                temp_static = static;
                temp_static(i) = temp_static(i)+h;
                f2 = eval_path_constraints(structure,x_in,x_0,x_f,t_0,t_f,temp_static,structure.uniform_els,0);
                Jstatic(i,:) = (f2-c)/h;
                
            end    

            
        else
            
            Jstatic = [];
            
        end
        
        % derivative wrt t0
        
        if any(structure.t0_vars)
            
            f2 = eval_path_constraints(structure,x_in,x_0,x_f,t_0+h,t_f,static,structure.uniform_els,0);
            dfdt0 = (f2-c)'/h;
            
        else
            
            dfdt0 = [];
            
        end
        
        % derivative wrt tf
        
        if any(structure.tf_vars)
            
            f2 = eval_path_constraints(structure,x_in,x_0,x_f,t_0,t_f+h,static,structure.uniform_els,0);
            dfdtf = (f2-c)'/h;
            
        else
            
            dfdtf = [];
            
        end
        
         % derivative wrt free x_0
        
        if n_x0>0
            
            dfdx0 = zeros(n_x0,size(Jc,2));
            
            for i =1:n_x0
                
                x0_temp = x_0;
                x0_temp(ids(i)) = x0_temp(ids(i))+h;
                f2 = eval_path_constraints(structure,x_in,x0_temp,x_f,t_0,t_f,static,structure.uniform_els,0);
                dfdx0(i,:) = (f2-c)'/h;
                
            end
            
        else
            
            dfdx0 = [];
            
        end

        % derivative wrt tf moved at the end, like the variable itself
        
        Jc = [Jstatic;dfdt0;dfdx0;Jc;dfdtf];
                
    end
    
end

% Boundary and integral inequality constraints

if ~isempty(structure.h)    
    
    [c2,Jc2] = eval_cost_functions2(structure.h,structure.wh,x_in,x_0,x_f,[t_0 t_f],static,structure.uniform_els,structure,jacflag,structure.dhu0,structure.dhxf,structure.dhuf,structure.dhxi,structure.dhui);
    c = [c;c2]; %add constraints
    
    % derivatives wrt generic static variables
    
    if jacflag
        
        if nstatic>0
            
            Jstatic = zeros(nstatic,size(Jc2,2));
            
            for i=1:nstatic
                
                temp_static = static;
                temp_static(i) = temp_static(i)+h;
                f2 = eval_cost_functions2(structure.h,structure.wh,x_in,x_0,x_f,[t_0 t_f],temp_static,structure.uniform_els,structure,0,[],[],[],[],[]);
                Jstatic(i,:) = (f2-c2)/h;
                
            end
            
        else
            
           Jstatic = zeros(0,length(c2));
            
        end
        
        % derivative wrt t0
        
        if any(structure.t0_vars)
            
            f2 = eval_cost_functions2(structure.h,structure.wh,x_in,x_0,x_f,[t_0+h t_f],static,structure.uniform_els,structure,0,[],[],[],[]);
            dfdt0 = (f2-c2)'/h;
            
        else
            
            dfdt0 = zeros(0,length(c2));
            
        end
        
        % derivative wrt final time
        
        if any(structure.tf_vars)
            
            f2 = eval_cost_functions2(structure.h,structure.wh,x_in,x_0,x_f,[t_0 t_f+h],static,structure.uniform_els,structure,0,[],[],[],[]);
            dfdtf = (f2-c2)'/h;
            
        else
            
            dfdtf = zeros(0,length(c2));
            
        end
        
        % derivative wrt free x_0
        
        if n_x0>0
            
            dfdx0 = zeros(n_x0,size(Jc2,2));
            
            for i =1:n_x0
                
                x0_temp = x_0;
                x0_temp(ids(i)) = x0_temp(ids(i))+h;
                f2 = eval_cost_functions2(structure.h,structure.wh,x_in,x0_temp,x_f,[t_0 t_f],static,structure.uniform_els,structure,0,[],[],[],[]);
                dfdx0(i,:) = (f2-c2)'/h;
                
            end
            
        else
            
            dfdx0 = zeros(0,length(c2));
            
        end
        
        % derivative wrt tf moved at the end, like the variable itself        
        
        Jc2 = [Jstatic;dfdt0;dfdx0;Jc2;dfdtf];
        
        Jc = [Jc Jc2];
        
    end
    
end

%% Equality constraints

% Differential path constraints (Dynamics)

[ceq,Jceq] = eval_constraints(structure,x_in,x_0,x_f,t_0,t_f,static,structure.uniform_els,jacflag);

if jacflag
    
    Jceq = Jceq';
    
    % derivatives wrt generic static variables
       
    if nstatic>0
        
        Jstatic = zeros(nstatic,size(Jceq,2));
        
        for i=1:nstatic
            
            temp_static = static;
            temp_static(i) = temp_static(i)+h;
            f2 = eval_constraints(structure,x_in,x_0,x_f,t_0,t_f,temp_static,structure.uniform_els,0);
            Jstatic(i,:) = (f2-ceq)/h;
            
        end
        
    else
        
        Jstatic = [];
        
    end
    
    % derivative wrt initial time
    
    if any(structure.t0_vars)
        
        f2 = eval_constraints(structure,x_in,x_0,x_f,t_0+h,t_f,static,structure.uniform_els,0);
        
        dfdt0 = (f2-ceq)'/h;
        
    else
        
        dfdt0 = [];
        
    end
    
    % derivative wrt final time
    
    if any(structure.tf_vars)
        
        f2 = eval_constraints(structure,x_in,x_0,x_f,t_0,t_f+h,static,structure.uniform_els,0);
        
        dfdtf = (f2-ceq)'/h;
        
    else
        
        dfdtf = [];
        
    end
    
    % derivative wrt free x_0
    
    if n_x0>0
        
        dfdx0 = zeros(n_x0,size(Jceq,2));
        
        for i =1:n_x0
            
            x0_temp = x_0;
            x0_temp(ids(i)) = x0_temp(ids(i))+h;
            f2 = eval_constraints(structure,x_in,x0_temp,x_f,t_0,t_f,static,structure.uniform_els,0);
            dfdx0(i,:) = (f2-ceq)'/h;
            
        end
        
    else
        
        dfdx0 = [];
        
    end
    
    % derivative wrt tf moved to the end, like the variable
    
    Jceq = [Jstatic;dfdt0;dfdx0;Jceq;dfdtf];
    
end

% Algebraic path constraints

if ~isempty(structure.e)
    
    [ceq2,Jceq2] = eval_path_eq_constraints(structure,x_in,x_0,x_f,t_0,t_f,static,structure.uniform_els,jacflag);
    ceq = [ceq;ceq2];
        
    % derivatives wrt generic static variables
    
    if jacflag

        Jceq2 = Jceq2';        
        
        if nstatic>0
            
            Jstatic = zeros(nstatic,size(Jceq2,2));
            
            for i=1:nstatic
                
                temp_static = static;
                temp_static(i) = temp_static(i)+h;
                f2 = eval_path_eq_constraints(structure,x_in,x_0,x_f,t_0,t_f,temp_static,structure.uniform_els,0);
                Jstatic(i,:) = (f2-ceq2)/h;
                
            end
            
        else
            
            Jstatic = [];
            
        end
        
        % derivative wrt t0
        
        if any(structure.t0_vars)
            
            f2 = eval_path_eq_constraints(structure,x_in,x_0,x_f,t_0+h,t_f,static,structure.uniform_els,0);
            dfdt0 = (f2-ceq2)'/h;
            
        else
            
            dfdt0 = [];
            
        end
        
        % derivative wrt tf
        
        if any(structure.tf_vars)
            
            f2 = eval_path_eq_constraints(structure,x_in,x_0,x_f,t_0,t_f+h,static,structure.uniform_els,0);
            dfdtf = (f2-ceq2)'/h;
            
        else
            
            dfdtf = [];
            
        end
        
        % derivative wrt free x_0
        
        if n_x0>0
            
            dfdx0 = zeros(n_x0,size(Jceq2,2));
            
            for i =1:n_x0
                
                x0_temp = x_0;
                x0_temp(ids(i)) = x0_temp(ids(i))+h;
                f2 = eval_path_eq_constraints(structure,x_in,x0_temp,x_f,t_0,t_f,static,structure.uniform_els,0);
                dfdx0(i,:) = (f2-ceq2)'/h;
                
            end
            
        else
            
            dfdx0 = [];
            
        end
                
        % derivative wrt tf moved to the end, like the variable
        
        Jceq2 = [Jstatic;dfdt0;dfdx0;Jceq2;dfdtf];
        
        Jceq = [Jceq Jceq2];
        
    end
    
end

% Boundary and integral constraints

if ~isempty(structure.q)
    
    % Final condition path constraints and Integral path constraints
    
    [ceq2,Jceq2] = eval_cost_functions2(structure.q,structure.wq,x_in,x_0,x_f,[t_0 t_f],static,structure.uniform_els,structure,jacflag,structure.dqu0,structure.dqxf,structure.dquf,structure.dqxi,structure.dqui);
    ceq = [ceq;ceq2];
    
    % derivatives wrt generic static variables
    
    if jacflag
        
        if nstatic>0
            
            Jstatic = zeros(nstatic,size(Jceq2,2));
            
            for i=1:nstatic
                
                temp_static = static;
                temp_static(i) = temp_static(i)+h;
                f2 = eval_cost_functions2(structure.q,structure.wq,x_in,x_0,x_f,[t_0 t_f],temp_static,structure.uniform_els,structure,0,[],[],[],[]);
                Jstatic(i,:) = (f2-ceq2)'/h;
                
            end
            
        else
            
            Jstatic = [];
            
        end
        
        % derivative wrt initial time
        
        if any(structure.t0_vars)
            
            f2 = eval_cost_functions2(structure.q,structure.wq,x_in,x_0,x_f,[t_0+h t_f],static,structure.uniform_els,structure,0,[],[],[],[]);
            dfdt0 = (f2-ceq2)'/h;
            
        else
            
            dfdt0 = [];
            
        end
        
        % derivative wrt final time
        
        if any(structure.tf_vars)
            
            f2 = eval_cost_functions2(structure.q,structure.wq,x_in,x_0,x_f,[t_0 t_f+h],static,structure.uniform_els,structure,0,[],[],[],[]);
            dfdtf = (f2-ceq2)'/h;
            
        else
            
            dfdtf = [];
            
        end
        
        % derivative wrt free x_0
        
        if n_x0>0
            
            dfdx0 = zeros(n_x0,size(Jceq2,2));
            
            for i =1:n_x0
                
                x0_temp = x_0;
                x0_temp(ids(i)) = x0_temp(ids(i))+h;
                f2 = eval_cost_functions2(structure.q,structure.wq,x_in,x0_temp,x_f,[t_0 t_f],static,structure.uniform_els,structure,0,[],[],[],[]);
                dfdx0(i,:) = (f2-ceq2)'/h;
                
            end
            
        else
            
            dfdx0 = [];
            
        end
        
        % derivative wrt tf moved to the end, like the variable
        
        Jceq2 = [Jstatic;dfdt0;dfdx0;Jceq2;dfdtf];
        
        Jceq = [Jceq Jceq2];
        
    end
    
end

end