function [val,dval] = join_objectives (objectives,grad_objectives,jacflag,problem)

% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2017 University of Strathclyde and Authors ------
%--------------- e-mail: lorenzo.ricciardi@strath.ac.uk----------------
%-------------------- Author: Lorenzo A. Ricciardi --------------------
%
% Calls the user defined function that links the objectives of the various
% phases. Also computes the Jacobian of this link function wrt the objectives 

val = problem.g(objectives);

%% Compute Jacobian of join function

if jacflag
    
    h = 1e-9;
    
    % get number of objective functions in total (summing all objs of all
    % phases) and create padding matrix with jacobians of objectives wrt states
    % and controls
    
    nobjs = 0;
    nvars = 0;
    
    for i = 1:length(objectives)
        
        nobjs = nobjs+length(objectives{i});
        nvars = nvars+size(grad_objectives{i},1);
        
    end
    
    padmatr = zeros(nvars,nobjs);
    onevect = ones(nvars,1);
    
    startx = 1;
    
    for i = 1:problem.num_phases%length(objectives)
        
        endx = startx+size(grad_objectives{i},2)-1;
        position = logical(onevect.*(problem.phase_mask==i));
        padmatr(position,startx:endx) = grad_objectives{i};
        startx = endx+1;
        
    end
    
    % Jacobian of join function wrt objectives
    
    gg = zeros(nobjs,length(val));
    
    posy = 1;    
    
    for i = 1:problem.num_phases                       
        
        for j = 1:length(objectives{i})
            
            tmpobj = objectives;            % clone objectives cell array
            incr = tmpobj{i};               % get vector of objectives of phase i
            incr(j) = incr(j)+h;            % increment jth objective of phase i
            tmpobj{i} = incr;               % re insert modified objective vector into objectives cell array
            
            tmpval = problem.g(tmpobj);     % call join function with modified objectives cell
            
            gg(posy,:) = (tmpval-val)'/h;       % compute finite difference approx
            posy = posy+1;
            
        end
        
    end
    
    % Jacobian of join function wrt actual optimisation variables (multiply by
    % precomputed gradients of individual objectives wrt their vars, pad
    % matrices)
    % NEED TO CHECK FOR MULTIPLE OBJECTIVES!!!!
    
    dval = padmatr*gg;
    
else
    
    dval = [];
    
end



end