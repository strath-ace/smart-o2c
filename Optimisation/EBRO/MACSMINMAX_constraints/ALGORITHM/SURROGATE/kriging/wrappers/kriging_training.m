function [model] = kriging_training(x, y, surrogate)

n = size(x,2);
if strcmp(func2str(surrogate.corrfun),'correxpg')
    n = n+1;
end

theta = repmat(1.0,1,n);
lob = repmat(1e-1,1,n);
upb = repmat(20,1,n);
model = dacefit(x, y, surrogate.regrfun, surrogate.corrfun, theta, lob, upb);
% model = dacefit(x, y, surrogate.regrfun, surrogate.corrfun, theta);

return