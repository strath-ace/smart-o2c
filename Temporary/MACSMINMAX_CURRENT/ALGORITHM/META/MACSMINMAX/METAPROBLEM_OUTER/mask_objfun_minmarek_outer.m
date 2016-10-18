function [masked] = mask_objfun_minmarek_outer(d,par_objfun)

[y,mse] = par_objfun.surrogate.predictor(d, par_objfun.surrogate.model);

masked = -par_objfun.surrogate.indicator(y,mse,par_objfun.ymin); % negative because we will always maximize PI or EI

return