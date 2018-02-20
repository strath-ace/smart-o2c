% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
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