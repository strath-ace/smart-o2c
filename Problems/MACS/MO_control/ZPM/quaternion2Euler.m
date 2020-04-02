function angles = quaternion2Euler(q)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%
angles = zeros(size(q,1),3);

for i = 1:size(q,1)
   
    q0 = q(i,1);
    q1 = q(i,2);
    q2 = q(i,3);
    q3 = q(i,4);
    
    angles(i,1) = atan2(2*(q0*q1+q2*q3),1-2*(q1*q1+q2*q2)); % phi
    angles(i,2) = asin(2*(q0*q2-q3*q1));                    % theta
    angles(i,3) = atan2(2*(q0*q3+q1*q2),1-2*(q2*q2+q3*q3)); % psi
    
end

end