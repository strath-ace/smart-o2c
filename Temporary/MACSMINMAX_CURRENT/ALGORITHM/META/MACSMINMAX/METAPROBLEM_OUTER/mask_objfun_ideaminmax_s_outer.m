function [masked] = mask_objfun_ideaminmax_s_outer(d,par_objfun)

n_obj = max(par_objfun.objectives);
if(n_obj>1)
    error('this needs to be checked because it does not work for MO!')
end

idx=ones(1,n_obj);
for obj = par_objfun.objectives
    for idx_u = 1:size(par_objfun.u_record{obj},1)
        [y,mse] = par_objfun.surrogate.predictor([d par_objfun.u_record{obj}(idx_u,:)], par_objfun.surrogate.model);
        if (idx_u == 1 || par_objfun.sign_inner*y > par_objfun.sign_inner*fmax)
            fmax = y;
            idx(obj) = idx_u;
        end
    end

end

u = par_objfun.u_record{1}(idx(1),:);
masked = -par_objfun.surrogate.indicator([d u],par_objfun.ymin,par_objfun.surrogate); % negative because we will always maximize PI or EI

return