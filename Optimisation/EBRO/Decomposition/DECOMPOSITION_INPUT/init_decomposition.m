% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function [in, problem] = init_decomposition()

% HERE YOU CAN DEFINE:

% 1) THE TIPE OF OUTPUT (Belief, Plausibility or both);
% 2) THE TIPE OF INPUT  (Run minmax and minmin  OR load them);
% 3) IF YOU WANT TO EVALUATE THE EXACT CURVE;
% 4) THE NUMBER OF SUB-SYSTEMS IN WHICH THE FUNCTION IS DECOMPOSED IN;
% 5) THE NUMBER OF SAMPLES: the bigger is the number, the more accurate is the decomposition;
% 6) THE DIMENTIONS OF COUPLED AND UNCOUPLED VECTORS;
% 7) LOWER AND UPPER BOUNDS OF DESIGN VECTOR (d) AND UNCERTAIN VECTOR (u), AND DIMENTION OF d;
% 8) NUMBER OF OBJECT FUNCTIONS;
% 9) FUNCTION(s);
% 10)MAX NUMBER OF EVALUATIONS FOR THE OPTIMIZER. 



%% 1) choose the output
% in.output = 0 --> only Belief (from min-max)
% in.output = 1 --> only Plausibility (from min-min)
% in.output = 2 --> Belief and Plausibility

in.output = 2;

%% 2) choose the input
% in.input = 0 --> do minmax and minmin
% in.input = 1 --> load d, u_min, u_max
% in.input = 2 --> load d, run max and min

in.input = 2;

%% 3) do exact Belief
% in.exact_curve(s) = 0 --> no
% in.exact_curve(s) = 1 --> yes

in.exact_curves = 1;

%% 4)number of sub-functions decomposition

num_functions = 2;  % number of sub-functions in which the problem is decomposed 

%% 5) number of samples
num_samples = 2;    % number of samples for each Belief and Plausibility curve of coupled vector





in.num_functions = num_functions;                                          
for i = 1:in.num_functions/2*(in.num_functions-1)
    in.num_samples{i}=num_samples;
end

%% EPISTEMIC VECTOR u (focal elements)

in.dim_u = [2 2 2]; %  6) dimentions of the vectors u_1, u_2, u_12, where:
%                         u_1 contains all the variables that influence only the sub-function f_1;
%                         u_2 contains all the variables that influence only the sub-function f_2;
%                         u_12 contains all the variables that influence both f_1 and f_2;

% 7) bounds of uncertain vector u;  

in.lb_u{1} = {[-5,-3,1]; [-5,-3,1]; [-5,-3,1]; [-5,-3,1]; [-5,-3,1]; [-5,-3,1] };   %  in.lower_bound_u{objective_function}
in.ub_u{1} = {[-1,0,2]; [-1,0,2]; [-1,0,2]; [-1,0,2]; [-1,0,2]; [-1,0,2] };  

in.bpa{1} = {[0.3 .3 0.4]; [0.3 .3 0.4]; [0.3 .3 0.4]; [0.3 0.3 0.4]; [0.3 0.3 0.4]; [0.3 0.3 0.4] };


%% DESIGN VECTOR d

dim_d = 2;           % 7) dimentions of the vector d, lower and upper bounds;

in.dim_d = dim_d;

in.lb_d = [-5;-5];
in.ub_d = [5;5];



%% FUNCTION 

% type of problem
problem.sign_inner = 1;               % -1 will run minmin(NOT CHANGE)

% objectives
problem.n_obj = 1;                    % 8) number of objective functions;
problem.objfun ={@function_F};        % 9) objective function(s); 
problem.par_objfun = {struct};

% design variables
                                                             
problem.dim_d = in.dim_d;                                           
problem.lb_d = in.lb_d;                                             
problem.ub_d = in.ub_d;                                             

% uncertain variables

problem.dim_u_i = in.dim_u;                                         
dim_u = sum(problem.dim_u_i);

problem.dim_u = dim_u;

for n =1:problem.n_obj

  
problem.lb_u{n} = in.lb_u{n}; 
problem.ub_u{n} = in.ub_u{n}; 


end

% 10) max number of optimizations

problem.maxnfeval = 50*5.0e3;

end

    