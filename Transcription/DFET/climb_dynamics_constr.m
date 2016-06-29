function [c,ceq,Jc,Jceq] = climb_dynamics_constr (f,g,h,wcis,structure,x_in,x_0,x_f,t_0,dfx,dfu)

% min time formulation
t_f = x_in(1);
x = x_in(2:end);

% max vel formulation
%t_f = 200;
%x = x_in;

%% Inequality constraints

c  = [];
Jc = [];

if ~isempty(g)
    
    % Normal path constraints
    
    [c,Jc] = eval_path_constraints(g,structure,x,x_0,x_f,t_0,t_f,structure.uniform_els,1,1,[],[]);
    
    Jc = Jc';
    
    % min time formulation needs derivative wrt final time, which is the
    % first optimisation variable
    
    f2 = eval_path_constraints(g,structure,x,x_0,x_f,t_0,t_f+1e-6,structure.uniform_els,1,0,[],[]);
    dfdt = (f2-c)'/1e-6;
    Jc = [dfdt;Jc];
    
end

if ~isempty(h)
    
    % Final condition path constraints and Integral path constraints
    
    [c2,Jc2] = eval_cost_functions2(h,wcis,x,x_f,[t_0 t_f],structure.uniform_els,structure,1,[],[]);
    
    Jc2 = Jc2';
    
    % min time formulation needs derivative wrt final time, which is the
    % first optimisation variable
    
    f2 = eval_cost_functions2(h,wcis,x,x_f,[t_0 t_f+1e-6],structure.uniform_els,structure,0,[],[]);
    dfdt = (f2-c2)'/1e-6;
    Jc2 = [dfdt;Jc2];
    
    c = [c;c2];
    Jc = [Jc Jc2];
    
end

%% Equality constraints

% Dynamics

[ceq,Jceq] = eval_constraints(f,structure,x,x_0,x_f,t_0,t_f,structure.uniform_els,1,dfx,dfu);

Jceq = Jceq';

% min time formulation needs derivative wrt final time, which is the
% first optimisation variablef2 = eval_constraints(f,structure,x,x_0,x_f,t_0,t_f+1e-6,structure.uniform_els,0,[],[]);
f2 = eval_constraints(f,structure,x,x_0,x_f,t_0,t_f+1e-6,structure.uniform_els,0,dfx,dfu);

dfdt = (f2-ceq)'/1e-6;
Jceq = [dfdt;Jceq];

end