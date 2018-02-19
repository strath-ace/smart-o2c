% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [surrogate] = surrogate_add(x_doe,f_doe,surrogate) % flag == 2

%function [surrogate] = surrogate_add_G(x_doe_aux,f_doe_aux,num_punti_doe,dim_u,dim_d,ZETA)    %% flag == 2

%         for i= 1:dim_u
%             num_intervalli(i) = map_info.n_int{i};
%             %num_extrem_division(i) = map_info.n_int{i}+1;                 
%         end

num_intervalli = surrogate.map_info.n_int;      
dim_d = surrogate.dim_d;
dim_u = surrogate.dim_u;

num_FE = size(surrogate.model,2);
surrogate.num_FE = num_FE;
d_doe = x_doe(:,1:dim_d);
u_doe = x_doe(:,dim_d+1:dim_d+dim_u); 
num_punti_doe = size(u_doe,1);
    
%% find FOCAL ELEMENT that contains the DOE's points

    for k = 1:num_punti_doe                          % line of matrix 'u_doe' = a single new point

        
        if size(u_doe,2) == 0 || num_FE == 0         % stop for-loop if I haven't FE or u(i). 
            break
        end
           

        n_FE = surrogate.find(x_doe(k,:),surrogate);

        % surrogate.x_doe(posizione*num_punti_doe-(num_punti_doe-1)+k-1,:) =  [d_doe(k,:) u_doe(k,:)];  % non necessario
        surrogate.changed_FE(n_FE) = surrogate.changed_FE(n_FE) + 1;

%% cell 
    
        
        surrogate.x_doe{n_FE} = [surrogate.x_doe{n_FE} ; d_doe(k,:) u_doe(k,:)];
        surrogate.f_doe{n_FE} = [surrogate.f_doe{n_FE} ; f_doe(k)];
    end
    
    
end