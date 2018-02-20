% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [n_FE] = surrogate_find_one(x_doe,surrogate)

% num_intervalli = surrogate.map_info.n_int;      
% dim_d = surrogate.dim_d;
% dim_u = surrogate.dim_u;

% d_doe = x_doe(1,1:dim_d);
% u_doe = x_doe(1,dim_d+1:dim_d+dim_u); 
    
% %% find FOCAL ELEMENT that contains the DOE's points
                 
% for i = 1:dim_u                              % number of the column of ZETA = u(i)
%     for j = 1:num_intervalli{i}              % da cmbiare sotto
%         if u_doe(1,i) <=  surrogate.map_info.interval_indicator{i}(1,j)  %map_info.interval_indicator{i}(j) 
%             vettore(i)=j;                    % j is the interval in u(i) where is u_doe(i)
%             break
%         end
%     end
% end  

% n_FE = 0;                                     % number(i) of focal element
% num_divisioni_vecchie = 1;

% for i = 1:dim_u-1 

%     num_divisioni_vecchie = num_divisioni_vecchie*num_intervalli{i};

%     n_FE = n_FE + (vettore(i+1)-1)*num_divisioni_vecchie;     % indici: matrice --> vettore num_divisioni_vecchie(i)

% end
% n_FE = n_FE + vettore(1); 
    

n_FE =1;

end