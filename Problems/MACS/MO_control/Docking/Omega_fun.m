function Omega = Omega_fun(w1,w2,w3)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%


Omega = [  0  w3 -w2 w1;
         -w3   0  w1 w2;
          w2 -w1   0 w3;
         -w1 -w2 -w3  0];

end