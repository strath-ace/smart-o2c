function [PC,DPC,elem_basis,basis_eval] = make_basis2(order,num_elems,num_eqs,x_temp,int_nodes,distrib)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%
% Generates polynomial basis: Bernstein basis or a Lagrange-interpolated
% basis on Legendre/Lobatto/Radau nodes

% As these are polynomials, only the coefficients are needed, and
% derivation/integration can be performed exactly :)

% Coefficients of polynomials of state
PC = zeros(order+1,num_elems*(order+1),num_eqs);
DPC = zeros(order,num_elems*(order+1),num_eqs);

% this is a selection vector, which will permute to build the base
% eg [1 0 0] -> [0 1 0] -> [0 0 1]

sel_vect = zeros(1,order+1);
sel_vect(1) = 1;
elem_basis = cell(num_elems,1);
basis_eval = cell(num_elems,length(int_nodes));

for i = 1:num_elems
    
    for k = 1:num_eqs
        
        for j =1:order+1
            
            if strcmp(distrib,'Legendre')||strcmp(distrib,'Lobatto')
                
                % this is just a simple method to get the coefficients, a more
                % thorough implementation could give the possibility of using
                % polynomials of order >20
                
                y_temp = circshift(sel_vect,[1,j-1]);
                PC(:,j+(i-1)*(order+1),k) = polyfit(x_temp,y_temp,order)';
                
                %cleanup of noise
                PC(abs(PC(:,j+(i-1)*(order+1),k))<1e-6,j+(i-1)*(order+1),k) = 0;
                
                DPC(:,j+(i-1)*(order+1),k) = polyder(PC(:,j+(i-1)*(order+1),k))';
                
            else
                
                % Bernstein
                
                PC(:,j+(i-1)*(order+1),k) = Bernstein_basis_coefficients(j,order)';
                
                % Coefficients of Bernstein Polynomials are integer :D
                % No noiste to clean
                
                DPC(:,j+(i-1)*(order+1),k) = polyder(PC(:,j+(i-1)*(order+1),k))';

                
            end
            
        end
        
    end
    
    elem_basis{i} = coeffs_to_horner(PC(:,1+(i-1)*(order+1):order+1+(i-1)*(order+1),:));
    
    for q = 1:length(int_nodes)
        
        % storing as sparse matrix, with trimmed near zero elements
        
        temp = elem_basis{i}(int_nodes(q));
        temp(abs(temp)<1e-9) = 0;
        temp = sparse(temp);
        
        basis_eval{i,q} = temp;
        
    end
    
end

end