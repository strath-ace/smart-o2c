% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%-----------Copyright (C) 2016 University of Strathclyde-------------
%


function [els, el_nodes, nodes] = make_elements(t_n,order,type)

num_elems = length(t_n)-1;

els = zeros(num_elems,2);
el_nodes = zeros(num_elems,order+1);

if strcmp(type,'Cheby')

tmpnodes = cos(fliplr((pi/(2*(order+1))):(pi/(order+1)):pi));

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