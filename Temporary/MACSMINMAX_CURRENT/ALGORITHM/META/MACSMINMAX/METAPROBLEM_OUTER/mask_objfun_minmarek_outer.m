function [masked] = mask_objfun_minmarek_outer(d,par_objfun)

masked = -par_objfun.surrogate.indicator(d,par_objfun.ymin,par_objfun.surrogate); % negative because we will always maximize PI or EI

return