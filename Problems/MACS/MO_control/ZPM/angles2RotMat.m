function A = angles2RotMat (angles)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%

    D = [cos(angles(1)) sin(angles(1)) 0; -sin(angles(1)) cos(angles(1)) 0; 0 0 1];
    C = [1 0 0; 0 cos(angles(2)) sin(angles(2)); 0 -sin(angles(2)) cos(angles(2))];
    B = [cos(angles(3)) sin(angles(3)) 0; -sin(angles(3)) cos(angles(3)) 0; 0 0 1];
    A = B*C*D;
    
end