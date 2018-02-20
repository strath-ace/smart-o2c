function h = coeffs_to_horner (coeffs)
<<<<<<< HEAD:Transcription/DFET/coeffs_to_horner.m

=======
>>>>>>> 5b7361d93c9119cf1d2e9e6c885bed93f924d71b:Transcription/DFET/coeffs_to_horner.m
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
<<<<<<< HEAD:Transcription/DFET/coeffs_to_horner.m
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
=======
%-----------Copyright (C) 2016 University of Strathclyde-------------
%
%
>>>>>>> 5b7361d93c9119cf1d2e9e6c885bed93f924d71b:Transcription/DFET/coeffs_to_horner.m
%
% Returns a function handle h for the evaluation of a polynomial with
% coeffs coefficient. The function is written in Horner's form for higher
% computational performance and accuracy.
% If coeffs is a row vector, it returns a scalar function handle
% If it is a matrix, generates a vector function handle
% If it is a 3d matrix, generates a matrix function handle

s = '@(t) [';
SC = size(coeffs);
num_coeffs = SC(1);
num_funcs = SC(2);

if length(SC)==3
    
    num_eqs = SC(3);
    
else
    
    num_eqs = 1;
    
end

for k = 1:num_eqs
    
    % leading zeroes
           
    s = strcat(s,repmat('0,',1,(k-1)*num_funcs));
        
    % Horner's expression for the polynomial

    for i = 1:num_funcs        
                
        for j=0:num_coeffs-2
            
            s = strcat(s,num2str(coeffs(end-j,i,num_eqs),16),'+t.*(');
            
        end
        
        s = strcat(s,num2str(coeffs(1,i,num_eqs),16));
        
        for j =0:num_coeffs-2
            
            s = strcat(s,')');
            
        end
        
        s = strcat(s,',');
        
    end
    
    % trailing zeroes
    
    s = strcat(s,repmat('0,',1,(num_eqs-k)*num_funcs-1));
    
    % last item with line change
    
    if k<num_eqs
        
        s = strcat(s,'0;');
        
    else
       
        s= strcat(s,';');
        
    end
    
end

s = strcat(s,']');

h = str2func(s);

end
