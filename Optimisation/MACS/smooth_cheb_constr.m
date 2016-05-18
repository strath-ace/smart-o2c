function [c,ceq]=smooth_cheb_constr(in,lambda,z,zstar,options,t_0,t_f,x_0,x_f)

alpha = in(1);
x = in(2:end);

% Equality constraints: dynamics. For now, no Jacobians are computed
% because I don't know if we can give Jacobians only of equality
% constraints and not of inequality constraints...

[~,ceq] = dynamics(options.oc.structure,x(options.oc.transcription_vars),x_0,x_f,t_0,t_f,0);

cc = eval_control_objectives(x(options.oc.transcription_vars),options.oc.structure,x_f,t_0,t_f);

% Re-scaled scalarisation
% BEWARE, THIS APPROACH FAILS FOR SINGLE OBJECTIVES (cc=z, 
% thus alpha must be equal to zero and thus the objective cannot change!)
% if any(lambda==1)
%     id = 1:length(lambda);
%     id = id(lambda==1);
%     c = (lambda(id).*(cc(id)-z(id))./(((zstar(id)-z(id))~=0)+((zstar(id)-z(id))==0))-alpha)';
% else
    c = (lambda.*(cc-z)./(((zstar-z).*((zstar-z)~=0))+((zstar-z)==0))-alpha)';
%end

return