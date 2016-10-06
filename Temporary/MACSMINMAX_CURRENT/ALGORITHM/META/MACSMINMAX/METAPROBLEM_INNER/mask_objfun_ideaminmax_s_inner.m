function [masked] = mask_objfun_ideaminmax_s_inner(u,par_objfun)

masked = -par_objfun.surrogate.indicator([par_objfun.d u],par_objfun.ymin,par_objfun.surrogate); % negative because we will always maximize PI or EI

return