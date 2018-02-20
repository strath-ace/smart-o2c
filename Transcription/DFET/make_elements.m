function [els, el_nodes, nodes] = make_elements(t_n,order,type)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%
% Generates nodes on the elements, to construct the associated bases.

% For Bernstein basis, which are not defined on nodes, for consistency with
% the rest of the code it simply generates equispaced nodes, but these are
% never used nor needed

num_elems = length(t_n)-1;

els = zeros(num_elems,2);
el_nodes = zeros(num_elems,order+1);

if strcmp(type,'Bernstein')     % Bernstein polynomials are actually defined in a different way

    tmpnodes = linspace(-1,1,order+1);

else
   
    if strcmp(type,'Legendre')
    
        tmpnodes = fliplr(lgwt(order+1,-1,1));
        
    else

        if strcmp(type,'Lobatto')
    
            tmpnodes = fliplr(lglnodes(order)');
            
        end
        
    end
    
end


for i = 1:num_elems
    
    els(i,:) = [t_n(i) t_n(i+1)];
    
    el_nodes(i,:) = (els(i,2)+els(i,1))/2 + (els(i,2)-els(i,1))/2* tmpnodes;
end

if any(isnan(el_nodes))
   
    el_nodes = els;
    
end

nodes = tmpnodes;

end