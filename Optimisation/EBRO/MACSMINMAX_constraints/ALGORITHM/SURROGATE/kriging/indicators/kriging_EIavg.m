% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [indicator] = kriging_EIavg(y, mse, ymin)

EIaug = kriging_EIaug(y, mse, ymin);
EIdom = kriging_EIdom(y, mse, ymin);

indicator = (max(EIaug,0)+max(EIdom,0))/2;

end
