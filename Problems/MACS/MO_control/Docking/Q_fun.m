function Q = Q_fun(q1,q2,q3,q4)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%


Q = [q1^2+q4^2-q2^2-q3^2 2*(q1*q2+q4*q3) 2*(q1*q3-q4*q2);
     2*(q1*q2-q4*q3) q2^2+q4^2-q1^2-q3^2 2*(q2*q3+q4*q1);
     2*(q1*q3+q4*q2) 2*(q2*q3-q4*q1) q3^2+q4^2-q1^2-q2^2];

end