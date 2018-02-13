function [masked] = mask_objfun_ideaminmax_s_outer(d,par_objfun)

n_obj = max(par_objfun.objectives);
y = nan(1,n_obj);
mse = nan(1,n_obj);
% idx=ones(1,n_obj);
for obj = par_objfun.objectives
    for idx_u = 1:size(par_objfun.u_record{obj},1)
        n_FE = par_objfun.surrogate{obj}.find([d par_objfun.u_record{obj}(idx_u,:)],par_objfun.surrogate{obj});
        [y_aux,mse_aux] = par_objfun.surrogate{obj}.predictor([d par_objfun.u_record{obj}(idx_u,:)], par_objfun.surrogate{obj}.model{n_FE});
        if (idx_u == 1 || par_objfun.sign_inner*y_aux > par_objfun.sign_inner*y(1,obj))
            y(1,obj) = y_aux;
            mse(1,obj) = mse_aux;
            % fmax = y_aux(1,obj);
            % idx(obj) = idx_u;
        end
    end

end

% x=[];
% for obj = 1:n_obj
%     x(obj,:) = [d, par_objfun.u_record{obj}(idx(obj),:)];
% end
% masked = -par_objfun.surrogate.indicator(x,par_objfun.ymin,par_objfun.surrogate);
masked = -par_objfun.indicator_d(y,mse,par_objfun.ymin);

return