function h = coeffs_to_handle (coeffs) 

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%
% Returns a function handle h for the evaluation of a polynomial with
% coeffs coefficient
% If coeffs is a row vector, it returns a scalar function handle
% If it is a matrix, generates a vector function handle

s = '@(t) [';
SC = size(coeffs);
nr = SC(1);
cl = SC(2);

for j = 1:nr

    for i=1:cl-1
    
        s = strcat(s,num2str(coeffs(j,i),16),'.*t.^',num2str(cl-i),'+');

    end
    
    s = strcat(s,num2str(coeffs(j,end),16),';');

    
end

s = strcat(s,']');


h = str2func(s);

end