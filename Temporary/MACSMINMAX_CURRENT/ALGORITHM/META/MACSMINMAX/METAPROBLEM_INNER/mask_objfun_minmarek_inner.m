function [masked] = mask_objfun_minmarek_inner(u,par_objfun)

masked = -par_objfun.surrogate.indicator(u,par_objfun.ymin,par_objfun.surrogate); % negative because we will always maximize PI or EI

return