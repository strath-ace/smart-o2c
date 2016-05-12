function [val,u_out,x_out] = fastest_climb_MACS_fmincon(x_in,f,x_0,x_f,t_0,structure,dfx,dfu,lb,ub,fminconoptions,varargin)

t_f = x_in(1);      % t_f is an optimisation variable!!!

if ~isempty(varargin)
    
    if ~all(all(isinf(varargin{2})))
        
        u_old = varargin{1};
        u_old = u_old(2:end);
        x_temp = varargin{2};
        
        x_guess  = update_guess(u_old,x_temp,structure);
        
    else
        
        x_guess = make_first_guess(f,x_0,t_0,x_in(1),x_in(2:end)',structure);
        
    end
    
else
    
    x_guess = make_first_guess(f,x_0,t_0,x_in(1),x_in(2:end)',structure);
    
end

x_guess = [t_f;x_guess];

[x_sol,~,exitflag] = fmincon(@(x) 1,x_guess,[],[],[],[],lb,ub,@(x) dynamics2(f,structure,x,x_0,x_f,t_0,dfx,dfu),fminconoptions);

t_f = x_sol(1);

g = @(x,u,t) [t_f 0];
weights = [1 0];

[x_out,u,x_b] = extract_solution(x_sol(2:end),structure,x_f);
u_out = [t_f u(:)'];

if exitflag==1
    
    val = eval_cost_functions(g,weights,x_out,u,x_b,[t_0 t_f],structure,0,[],[]);
    
else
    
    val = inf;
    x_out = inf(size(x_out,1),size(x_out,2));
    
end

end

function [c,ceq,Jc,Jceq] = dynamics2(f,structure,x_in,x_0,x_f,t_0,dfx,dfu)

t_f = x_in(1);
x_guess = x_in(2:end);

c = [];
Jc = [];

[ceq,Jceq] = eval_constraints(f,structure,x_guess,x_0,x_f,t_0,t_f,1,dfx,dfu);

Jceq = Jceq';
f2 = eval_constraints(f,structure,x_guess,x_0,x_f,t_0,t_f+1e-6,0,[],[]);
dfdt = (f2-ceq)'/1e-6;
Jceq = [dfdt;Jceq];

end

% function [val] = objectives2 (x_in,t_0,x_f,structure)
% 
% t_f = x_in(1);      % t_f is an optimisation variable!!!
% x_sol = x_in(2:end);
% 
% g = @(x,u,t) [t_f 0];
% weights = [1 0];
% 
% [x,u,x_b] = extract_solution(x_sol,structure,x_f);
% 
% val = eval_cost_functions(g,weights,x,u,x_b,[t_0 t_f],structure,0,[],[]);
% 
% end