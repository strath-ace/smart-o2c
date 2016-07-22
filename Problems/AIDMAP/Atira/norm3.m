%% norm3: This function calculates the norm of a vector with 3 elements
%
%% Inputs:
% * vec         : The vector of which the norm is to be calculated
%
%% Outputs: 
% * norm        : The norm of the vector [real number]
%
%% Author(s): Juan Manuel Romero Martin (2014)
%  Email:  juan.romero-martin@strath.ac.uk

function norm = norm3(vec)

    % Compute the norm of a 3-dimensional vector
    norm = sqrt(vec(1)^2 + vec(2)^2 + vec(3)^2);

end

