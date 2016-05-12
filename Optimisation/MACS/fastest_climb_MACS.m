function [val,u_out,x_out] = fastest_climb_MACS(x_in,f,x_0,x_f,t_0,structure,tol_conv,maxits,dfx,dfu,lb,ub,varargin)

t_f = x_in(1);      % t_f is an optimisation variable!!!

%x_guess = make_first_guess(f,x_0,t_0,x_in(1),x_in(2:end)',structure);

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

[x_sol,resids]=solve_constraints(f,x_guess,x_0,x_f,t_0,t_f,structure,tol_conv,maxits,dfx,dfu,lb,ub);    % solution of underdetermined system

g = @(x,u,t) [t_f 0; 0 u(1)];
weights = [1 0; 0 1];

[x_out,u,x_b] = extract_solution(x_sol,structure,x_f);
u_out = [t_f u(:)'];
% loglog(1:length(resids),resids,1:length(resids),tol_conv*ones(1,length(resids)));
% 
% drawnow
if resids(end)<tol_conv
    
    val = eval_cost_functions(g,weights,x_out,u,x_b,[t_0 t_f],structure,0,[],[]);
    
else
    
    val = inf;
    x_out = inf(size(x_out,1),size(x_out,2));
    
end