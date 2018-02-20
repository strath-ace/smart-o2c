% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [num_FE] = number_Focal_Element_TEST(problem)

num_FE = 1;
for i = 1:problem.dim_u
    
    num_FE = num_FE*length(problem.lb_u{1}{i});
    
end

end