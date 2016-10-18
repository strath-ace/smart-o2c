function [masked] = mask_objfun_ideaminmax_s_inner(u,par_objfun)

[y_aux,mse_aux] = par_objfun.surrogate.predictor([par_objfun.d u], par_objfun.surrogate.model);
y = y_aux(1,par_objfun.objective);
mse = mse_aux(1,par_objfun.objective);

masked = -par_objfun.surrogate.indicator(y,mse,par_objfun.ymin); % negative because we will always maximize PI or EI

return