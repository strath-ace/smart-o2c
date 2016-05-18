function [x,resids,tvect]=solve_constraints(f,x,x_0,x_f,t_0,t_f,structure,tol,maxits,dfx,dfu,lb,ub)

% Newton's Method for Under-determined systems
% Author: Joseph Simonis
% Latest update: 03-01-06
%
% x_0 initial guess of solution
% structure contains strcture of the problem, necessary to evaluate F and J
% tol solution tolerance
% maxits maximum number of nonlinear iterations
%
%maxstep = 0.1;
[F,J]=eval_constraints(f,structure,x,x_0,x_f,t_0,t_f,1,dfx,dfu);
% [F2,J2]=eval_constraints(f,structure,x,x_0,x_f,t_0,t_f,1,[],[]);
% if max(max(abs(J2-J)))>1e-6
%     error('Check Jacobian...')
% end
%residual=norm(F);
residual=max(abs(F));
its=1;
resids(its,1)=residual;
%fprintf('\nIt.No.   ||F(u)||\n');
%fprintf('  %d   %e   \n', 0,residual);
t=1;
tvect(its,1) = t;
while(residual > tol && its<maxits)
    
    s=-pinv(full(J))*F;             % Solve the under-determined lin. sys. (maybe it could be reformulated, avoinding the direct and probably expensive computation of the pseudoinverse)
    %     norms = max(s);
    %     s = s*(norms<maxstep)+s/norms*maxstep*(norms>=maxstep); % allow a maximum of maxstep (my personal improvement)
    
    % bound constrain
%     if ~isempty(ub)
%         
%         [t,s] = alpha_clip2(x,1,s,lb,ub);
%         
%         % if t is too small is better to just use hard clipping
%         t = max(t,1e-6);
%         
%     end
    
    x=x+s*t;                          % Update x
    % hard clipping afterwards
    x(x>ub)=ub(x>ub);
    x(x<lb)=lb(x<lb);        
    
    if any(x>ub) || any(x<lb)
        
        error('Couldn''t enforce box constraints');
        
    end
    
    %linnorm=norm(J*s+F);
    [F,J]=eval_constraints(f,structure,x,x_0,x_f,t_0,t_f,1,dfx,dfu);
    %     [F2,J2]=eval_constraints(f,structure,x,x_0,x_f,t_0,t_f,1,[],[]);
    %     if max(max(abs(J2-J)))>1e-6
    %         error('Check Jacobian...')
    %     end
    %residual=norm(F);
    residual=max(abs(F));
    %fprintf('  %d   %e\n', its,residual);
    its=its+1;
    resids(its,1)=residual;
    tvect(its,1)=t;

end
