
function [f] = mask(u,d,objfun,par_objfun)

f = -objfun(d,u,par_objfun);

return