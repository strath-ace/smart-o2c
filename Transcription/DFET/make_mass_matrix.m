function M = make_mass_matrix(num_elems,state_order,test_order,num_eqs,PC,DPC,pos,neg,DFET)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%
% Generates "mass" matrix, i.e. the matrix of the integrals of the products
% between the time derivative of the test functions and the basis functions
% of the states. Since all bases are polynomial and the coefficients are
% available, this operation can be performed exactly, no numerical
% integration is required

% TODO: improve readability by factoring out terms indexing the matrix
% TODO: make it faster... 4 nested loops scream for vectorisation!

M = sparse (num_elems*(test_order+1)*num_eqs,num_elems*(state_order+1)*num_eqs);

% This is a diagonal block matrix, and if the polynomials are defined from
% -1 to 1 the blocks depend only on the degree of the polynomial. Thus,
% they can be computed offline once and for all and stored in a separate
% read only file
% In case f(x) = k*x, the RHS vector becomes LHS MATRIX, with interesting
% consequences on M

if DFET == 0    % FET formulation, M contains integrals of phi*d(phi)/dt
    
    for i = 1:num_elems
        
        for q = 1:num_eqs
            
            for j=1:test_order+1
                
                for k=1:state_order+1
                    
                    % conv computes the coefficients of the product of two polys
                    prod = conv(PC(:,j+(i-1)*(state_order+1),q)',DPC(:,k+(i-1)*(state_order+1),q)');
                    
                    % this computes the coefficients of the indefinite integral of a poly
                    inte_coeffs = polyint(prod);
                    
                    % thie evaluates the integral and places it in the right
                    % position in the matrix
                    M(j+(q-1)*(state_order+1)+(state_order+1)*(i-1)*num_eqs,k+(q-1)*(state_order+1)+(state_order+1)*(i-1)*num_eqs) = inte_coeffs*(pos-neg);
                    
                end
                
            end
            
        end
        
    end
    
else            % DGFET formulation, M contains integrals of d(phi)/dt*phi
    
    for i = 1:num_elems
        
        for q = 1:num_eqs
            
            for j=1:test_order+1
                
                for k=1:state_order+1
                    
                    % conv computes the coefficients of the product of two polys
                    prod = conv(DPC(:,j+(i-1)*(test_order+1),q)',PC(:,k+(i-1)*(state_order+1),q)');
                    
                    % this computes the coefficients of the indefinite integral of a poly
                    inte_coeffs = polyint(prod);
                    
                    % thie evaluates the integral and places it in the right
                    % position in the matrix
                    M(j+(q-1)*(test_order+1)+(test_order+1)*(i-1)*num_eqs,k+(q-1)*(state_order+1)+(state_order+1)*(i-1)*num_eqs) = inte_coeffs*(pos-neg);
                    
                end
                
            end
            
        end
        
    end
    
end

%% Modification of M to enforce ICs and/or match elements

% TODO: improve readability by factoring out terms indexing the matrix

if DFET == 0    % FET formulation, strong matching imposed
    
    % Enforcing IC
    
    for q = 1:num_eqs
        
        M(1+(state_order+1)*(q-1),1+(state_order+1)*(q-1):state_order+1+(state_order+1)*(q-1)) =  (PC(:,1:state_order+1,q)'*neg(size(PC,1)+1*(state_order==0):end))';
        
    end
    
    % Match elements (end condition of elem. i-1 equals initial condition
    % of elem i). End condition in CONTINUOUS GALERKIN is the value of the
    % interpolating polynomials at the right node (1)
    
    for i = 2:num_elems
        
        for q = 1:num_eqs
            
            M(1+(q-1)*(state_order+1)+(state_order+1)*num_eqs*(i-1),1+(q-1)*(state_order+1)+(state_order+1)*num_eqs*(i-1):state_order+1+(q-1)*(state_order+1)+(state_order+1)*num_eqs*(i-1)) = (PC(:,1+(i-1)*(state_order+1):state_order+1+(i-1)*(state_order+1),q)'*neg(size(PC,1)+1*(state_order==0):end))';
            M(1+(q-1)*(state_order+1)+(state_order+1)*num_eqs*(i-1),1+(q-1)*(state_order+1)+(state_order+1)*num_eqs*(i-2):state_order+1+(q-1)*(state_order+1)+(state_order+1)*num_eqs*(i-2)) = -(PC(:,1+(i-1)*(state_order+1):state_order+1+(i-1)*(state_order+1),q)'*pos(size(PC,1)+1*(state_order==0):end))';
            
        end
        
    end
    
else           % DGFET formulation, weak matching imposed
    
    % Match elements (end condition of elem. i-1 equals initial condition
    % of elem i). End condition in BI-DISCONTINUOUS GALERKIN is NOT the
    % value of the interpolating polynomials at the right node, BUT the
    % value of the end node (an unknown to be computed). HOWEVER, these end
    % nodes "cancel out", so the matrix is really shorter than the computed
    % one... Basically, the first row of each element's submatrix is moved
    % one row up
    
    if num_elems>1
        
        M2 = sparse(((test_order+1)+(test_order)*(num_elems-1))*num_eqs,size(M,2));
        
        % copy first element as is               
        
        startx = 1;
        starty = 1;
        startx_aug = 1;
        starty_aug = 1;
        
        endx = startx + (state_order+1)*num_eqs-1;
        endy = starty + (test_order+1)*num_eqs-1;
        endx_aug = startx_aug + (state_order+1)*num_eqs-1;
        endy_aug = starty_aug + (test_order+1)*num_eqs-1;
        
        M2(starty_aug:endy_aug,startx_aug:endx_aug) = M(starty:endy,startx:endx);
        
        startx = endx+1;
        starty = endy+1;
        startx_aug = endx_aug+1;
        starty_aug = endy_aug+1;  % slip "up" by one
        
        pos_prev = 0;
        
        for i = 2:num_elems   %last element is different because it has final terms, so is treated separately
                        
            for q = 1:num_eqs
                
                pos_prev = pos_prev+(test_order+1)*((i-1)==1)+(test_order)*((i-1)>1);

                
                endx = startx+state_order;
                endy = starty+test_order;
                endx_aug = startx_aug+test_order-1;
                endy_aug = starty_aug+test_order-1;           
                    
                % copy this element submatrix into final matrix
                M2(starty_aug:endy_aug,startx_aug:endx_aug) = M(starty+1:endy,startx:endx);
                
                %impose matching with previous element
                M2(pos_prev,startx:endx_aug) = M(starty,startx:endx);
                
                startx = endx+1;
                starty = endy+1;
                startx_aug = endx_aug + 1;
                starty_aug = endy_aug + 1;
                
                
                
            end
            
        end
        
        M = M2;
        
    end
    
end

end