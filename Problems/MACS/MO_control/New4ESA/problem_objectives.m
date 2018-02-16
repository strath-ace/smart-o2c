function [val] = problem_objectives(objs)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%

% objs{1} : integral of uy^2 (phase 1)
% objs{2} : integral of ux^2 (phase 2)

o1 = objs{1};
o2 = objs{3};

%alpha = 0;

val =[o1(1);o2(1)];    

%val = objs{1};

end
