<<<<<<< HEAD:Transcription/DFET/extract_solution.m
function [x,u,xb] = extract_solution(x_sol,structure,xf)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%
% Reshapes solution vector into a state vector, a control vector, 
% and a boundary values vector 

x = zeros((structure.state_order+1)*structure.num_eqs,structure.num_elems);    %coefficients of the polynomials, elementwise
u = zeros((structure.control_order+1)*structure.num_controls,structure.num_elems);     %coefficients of the controls, elementwise

startx = 1;

for i = 1:structure.num_elems

    endx = startx+(structure.state_order+1)*structure.num_eqs-1;
    x(:,i) = x_sol(startx:endx);

    startx = endx+1*(structure.num_controls>0);
    
    endx = startx+(structure.control_order+1)*structure.num_controls-1*(structure.num_controls>0);
    
    % ugly, but apparently no empty to empty assignment can be done
    
    if structure.num_controls>0
    
        u(:,i) = x_sol(startx:endx);
        
    end
    
    startx = endx+1;
end

if structure.DFET==1
    
    count = 0;
    
    xb = zeros(size(xf));
    
    for i = 1:structure.num_eqs

        count = count+1*(structure.imposed_final_states(i)==0);        
        xb(i) = x_sol(startx-1+count)*(structure.imposed_final_states(i)==0)+xf(i)*(structure.imposed_final_states(i)==1);
        
    end
    
else
   
    xb = [];
    
end

end
=======
function [x,u,xb] = extract_solution(x_sol,structure,xf)
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%-----------Copyright (C) 2016 University of Strathclyde-------------
%
%
%

x = zeros((structure.state_order+1)*structure.num_eqs,structure.num_elems);    %coefficients of the polynomials, elementwise
u = zeros((structure.control_order+1)*structure.num_controls,structure.num_elems);     %coefficients of the controls, elementwise

startx = 1;

for i = 1:structure.num_elems

    endx = startx+(structure.state_order+1)*structure.num_eqs-1;
    x(:,i) = x_sol(startx:endx);

    startx = endx+1*(structure.num_controls>0);
    endx = startx+(structure.control_order+1)*structure.num_controls-1*(structure.num_controls>0);
    
    % ugly, but apparently no empty to empty assignment can be done
    
    if structure.num_controls>0
    
        u(:,i) = x_sol(startx:endx);
        
    end
    
    startx = endx+1;
end

if structure.DFET==1
    
    count = 0;
    
    for i = 1:structure.num_eqs

        count = count+1*(structure.imposed_final_states(i)==0);        
        xb(i) = x_sol(startx-1+count)*(structure.imposed_final_states(i)==0)+xf(i)*(structure.imposed_final_states(i)==1);
        
    end
    
else
   
    xb = [];
    
end

end
>>>>>>> 5b7361d93c9119cf1d2e9e6c885bed93f924d71b:Transcription/DFET/extract_solution.m
