% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [map_info] = get_map_info(lb_input,ub_input)
% assumes that the search space is the cartesian product of the intervals provided
% needs to be reimplemented for correlated variables to "prune out products"
% (e.g. include [1,2]x[1,2]U[3,4]x[3,4] but exclude [1,2]x[3,4])

dim = length(lb_input);
if (dim ~= length(ub_input))
    error('supplied incoherent bpa structure');
end

%outputs
map_info.scale = zeros(1,dim);
map_info.n_int = cell(dim,1);
map_info.lb = cell(dim,1);
map_info.interval_indicator = cell(dim,1);

for d = 1:dim
    lbi_input = lb_input{d};
    ubi_input = ub_input{d};
    n_int = length(lbi_input);
    if (n_int ~= length(ubi_input) || any(lbi_input > ubi_input))
        error('supplied incoherent bpa structure');
    end
    [lbi_input, idx] = sort(lbi_input);
    ubi_input = ubi_input(idx);

    lbi_output = lbi_input(1);
    ubi_output = ubi_input(1);
    for int = 2:n_int
        if (lbi_input(int) <= ubi_output(end))          % connected
            if (ubi_input(int) > ubi_output (end))      % but not completely contained
                ubi_output(end) = ubi_input(int);       % extend the previous output interval
                                                        % else (if completely contained) do nothing
            end
        else                                            % disconnected
            lbi_output = [lbi_output, lbi_input(int)];  % add another output interval
            ubi_output = [ubi_output, ubi_input(int)];
        end
    end

    interval_lengths = ubi_output - lbi_output;
    scale = sum(interval_lengths);
    interval_indicator = cumsum(interval_lengths)/scale;

    map_info.scale(d) = scale;
    map_info.n_int{d} = length(interval_indicator);
    map_info.interval_indicator{d} = interval_indicator;
    map_info.lb{d} = lbi_output;
end

return
