function out = C3 (q)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%

out = [2*( q(2)*q(4) - q(1)*q(3) );
       2*( q(3)*q(4) + q(1)+q(2) );
       1 - 2*(q(2)^2+q(3)^2)     ];

end