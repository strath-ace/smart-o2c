function [masked] = mask_objfun_mo_sa_outer(d,par_objfun)

[y,mse] = par_objfun.surrogate.predictor(d, par_objfun.surrogate.model);

masked = y;

return