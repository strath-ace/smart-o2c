function [X] = angles2Vectors(angles)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%

X = zeros(size(angles,1),9);

for i = 1:size(angles,1)
   
    M = angles2RotMat(angles(i,:));
    X(i,1:3) = (M*[1;0;0])';
    X(i,4:6) = (M*[0;1;0])';
    X(i,7:9) = (M*[0;0;1])';
    
end

end