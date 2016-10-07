function [masked] = mask_objfun_ideaminmax_s_outer(d,par_objfun)

n_obj = max(par_objfun.objectives);
idx=ones(1,n_obj);
for obj = par_objfun.objectives
    for idx_u = 1:size(par_objfun.u_record{obj},1)
        [y,mse] = par_objfun.surrogate.predictor([d par_objfun.u_record{obj}(idx_u,:)], par_objfun.surrogate.model);
        if (idx_u == 1 || par_objfun.sign_inner*y(1,obj) > par_objfun.sign_inner*fmax)
            fmax = y(1,obj);
            idx(obj) = idx_u;
        end
    end

end

x=[];
for obj = 1:n_obj
    x(obj,:) = [d, par_objfun.u_record{obj}(idx(obj),:)];
end
masked = -par_objfun.surrogate.indicator(x,par_objfun.ymin,par_objfun.surrogate);

% if n_obj == 1
% 	u = par_objfun.u_record{1}(idx(1),:);
% 	masked = -par_objfun.surrogate.indicator([d u],par_objfun.ymin,par_objfun.surrogate); % negative because we will always maximize PI or EI
% else
% 	error('this needs to be checked because it is not straightforward for MO!')
% 	% I think that a solution is to compute EIaug or EIdom with the y and mse values of a point (d,u) for each objective, where
% 	% the d would be the same. Conceptually this is like computing the EI for the minmax front that is kept in ymin, i.e. we are driving
% 	% the optimisation with an EI(d) even though we are using a single surrogate. That could work...
% end
return