% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%-----------Copyright (C) 2016 University of Strathclyde-------------
%


%% norm3: This function calculates the norm of a vector with 3 elements
% 
%% Inputs:
% * vec         : The vector of which the norm is to be calculated
% 
%% Outputs: 
% * norm        : The norm of the vector [real number]
% 
%% Author(s): Juan Manuel Romero Martin (2014)
% Email:  juan.romero-martin@strath.ac.uk

function norm = norm3(vec)

    % Compute the norm of a 3-dimensional vector
    norm = sqrt(vec(1)^2 + vec(2)^2 + vec(3)^2);

end

