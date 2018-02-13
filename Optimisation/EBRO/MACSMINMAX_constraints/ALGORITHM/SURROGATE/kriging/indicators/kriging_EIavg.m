function [indicator] = kriging_EIavg(y, mse, ymin)

EIaug = kriging_EIaug(y, mse, ymin);
EIdom = kriging_EIdom(y, mse, ymin);

indicator = (max(EIaug,0)+max(EIdom,0))/2;

end
