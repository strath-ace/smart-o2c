% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [lower_bound, upper_bound] =  getimit()
lower_bound = [-1000	3	0	0	100	100	30	400	800	0.01	0.01	...
    0.01	0.01	0.01	1.05	1.05	1.15	1.7	-pi -pi -pi -pi];

upper_bound = [0	5	1	1	400	500	300	1600	2200	0.9	0.9	0.9	...
    0.9	0.9	6	6	6.5	291	pi	pi	pi	pi];


