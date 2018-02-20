% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [Sample] = Sampling_Belief_Plausibility(Partial_curve, decomposition, in)

%
% SAMPLE each partial curve with the input-number of samples and save
% information about the vector u, boudaries of the Focal Elements, valoue
% of F and position i the partial curves.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_coupled_vector_tot = length(in.dim_u_i(in.num_functions+1:end));

for num_coupled_vector = 1:num_coupled_vector_tot
    
    
    if in.dim_u_i(in.num_functions + num_coupled_vector)>0
        num_samples_tot = in.num_samples{num_coupled_vector};
        
        
%% BELIEF
        
        if in.output == 0 || in.output == 2
            
            
            [Bel_unique_values, index] = sort(Partial_curve{num_coupled_vector}.Belief_FE_belief_partial(:));
            
            
            if num_samples_tot > length(Bel_unique_values)
                
                num_samples_tot = length(Bel_unique_values);          
            end
            
            Bel_subset_num = floor((length(Bel_unique_values)-1)/(num_samples_tot-1));
            
            
            
            
            
            for i = 1:num_samples_tot

                N_Focal_Element_coupled_vector_k = length(decomposition{1,num_coupled_vector}.FocalElement);
                
                if i==1
                    Sample{num_coupled_vector}.Belief_f(num_samples_tot+1-i) =  1;                                    
                    fun(i) = Partial_curve{num_coupled_vector}.Belief_FE_function(index(end)); 
                    Sample{num_coupled_vector}.position(num_samples_tot+1-i) = N_Focal_Element_coupled_vector_k;
                else
                    
                    

%                     Bel_subset = Bel_unique_values(end-(i-1)*Bel_subset_num : end-(i-2)*Bel_subset_num -1);
%                     [Sample{num_coupled_vector}.Belief_f(num_samples_tot+1-i), idx] = datasample(Bel_subset, 1);
%                     
%                     fun(i) = Partial_curve{num_coupled_vector}.Belief_FE_function(index(end-(i-1)*Bel_subset_num + idx -1));
%                     
%                     Sample{num_coupled_vector}.position(num_samples_tot+1-i) = index(end-(i-1)*Bel_subset_num + idx -1);% - par_treshold; par_treshold;% index(par_treshold);% 
% global par_treshold

% the following two lines make the samples from the right to the left and with step =1 over the focal elements.                    
                    Sample{num_coupled_vector}.Belief_f(num_samples_tot+1-i) = Partial_curve{num_coupled_vector}.Belief_FE_belief_partial(index(end-i+1)); %(index(par_treshold));%
                    fun(i) = Partial_curve{num_coupled_vector}.Belief_FE_function(index(end-i+1)); %(index(par_treshold));% sample in the right
                    Sample{num_coupled_vector}.position(num_samples_tot+1-i) = index(end-i+1);%index(par_treshold);%
                end
 



                
                for j = 1 : N_Focal_Element_coupled_vector_k
                    
                    % Belief
                    if abs(decomposition{1,num_coupled_vector}.FocalElement{j}.upper_f - fun(i)) == 0 
                        
                        Sample{num_coupled_vector}.FE{num_samples_tot+1-i} = decomposition{1,num_coupled_vector}.FocalElement{j};
                        
                    end
                    
                    
                end
                
            
            end
            
            Sample{num_coupled_vector}.f = fliplr(fun);
            
            Sample{num_coupled_vector}.num_coupled_vector = num_coupled_vector;

            
        end
% par_treshold = par_treshold +1;             
        
%% PLAUSIBILITY
        if in.output == 1 || in.output == 2
            
            
            
            [Plausibility_unique_values, index_Pl] = sort(Partial_curve{num_coupled_vector}.Plausibility_FE_plausibility_partial(:));
            
            if num_samples_tot >  length(Plausibility_unique_values)
                
                num_samples_tot = length(Plausibility_unique_values);          
            end
            
            Plausibility_subset_num = floor((length(Plausibility_unique_values)-1)/(num_samples_tot-1));
            
            
            
            for i = 1:num_samples_tot
                
                if i==1
                    
                    Sample{num_coupled_vector}.Plausibility_f(num_samples_tot+1-i) =  1;
                    fun_Pl(i) = Partial_curve{num_coupled_vector}.Plausibility_FE_function(index_Pl(end));
                else
                    
                    Plausibility_subset = Plausibility_unique_values(end-(i-1)*Plausibility_subset_num : end-(i-2)*Plausibility_subset_num -1);
%                     [Sample{num_coupled_vector}.Plausibility_f(num_samples_tot+1-i), idx_Pl] = datasample(Plausibility_subset, 1);
%                     
%                     fun_Pl(i) = Partial_curve{num_coupled_vector}.Plausibility_FE_function(index_Pl(end-(i-1)*Plausibility_subset_num + idx_Pl -1));
                    
% the following two lines make the samples from the right to the left and with step =1 over the focal elements.                    
                    Sample{num_coupled_vector}.Plausibility_f(num_samples_tot+1-i)=Partial_curve{num_coupled_vector}.Plausibility_FE_plausibility_partial(index_Pl(end-i+1));
                    fun_Pl(i) = Partial_curve{num_coupled_vector}.Plausibility_FE_function(index_Pl(end-i+1)); % sample in the right
                end
                
                for j=1:length(decomposition{1,num_coupled_vector}.FocalElement)
                    
                    if abs(decomposition{1,num_coupled_vector}.FocalElement{j}.downer_f - fun_Pl(i)) == 0 %< (minmax.f-minmin.f)/1000
                        
                        Sample{num_coupled_vector}.FE_Plausibility{num_samples_tot+1-i} = decomposition{1,num_coupled_vector}.FocalElement{j};
                        
                    end
                    
                end
                
                
            end
            
            Sample{num_coupled_vector}.f_Pl = fun_Pl;
            
            Sample{num_coupled_vector}.num_coupled_vector = num_coupled_vector;
            
        end
    
    end
    
end

end