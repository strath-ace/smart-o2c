% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [LIST] = plot_Belief_Plausibility_decomposition(decomposition_end, in)



fmax    = [];
bpa_max = [];
fmin    = [];
bpa_min = [];

fmax_add = [];
bpa_add  = [];

for i = 1:length(decomposition_end)
    if isempty(decomposition_end{i})==0
        
        if in.output == 0 || in.output == 2            
            % BELIEF           
            for k = 1 : length(decomposition_end{i}.step_two)     % vector of all the f-values of the FEs                
                fmax    = [fmax     decomposition_end{i}.step_two{k}.F]; %decomposition_end{i}.step_two{k}.F;
                bpa_max = [bpa_max  decomposition_end{i}.step_two{k}.bpa]; 
            end    
                %%
%                 fmax_add = [fmax_add decomposition_end{i}.FE_add.F];
%                 bpa_add  = [bpa_add  decomposition_end{i}.FE_add.bpa];
        end
        
        
        if in.output == 1 || in.output == 2            
            % PLAUSIBILITY            
            for k = 1 : length(decomposition_end{i}.step_two)     % vector of all the f-values of the FEs                
                fmin    = [fmin     decomposition_end{i}.step_two{k}.F_Plausibility];
                bpa_min = [bpa_min  decomposition_end{i}.step_two{k}.bpa_Plausibility];
            end           
        end
                
    end
end



%  Belief
if in.output == 0 || in.output == 2

% add the focal elements from an other threscold but with F<
    fmax    =  [fmax     fmax_add];
    bpa_max =  [bpa_max  bpa_add];
    
    %%
    [f_sorted_Bel, position_sorted_Bel] = sort(fmax);
    
    LIST.F_Bel = f_sorted_Bel;
    LIST.Bel   = cumsum(bpa_max(position_sorted_Bel));
    
    %%

end


% Plausibility
if in.output == 1 || in.output == 2
    
    [f_sorted_Pl, position_sorted_Pl] = sort(fmin);
    
    LIST.F_Pl = f_sorted_Pl;
    LIST.Pl   = cumsum(bpa_min(position_sorted_Pl));
end





end