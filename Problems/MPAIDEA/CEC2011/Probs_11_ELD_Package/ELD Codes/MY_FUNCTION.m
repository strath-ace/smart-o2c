<<<<<<< HEAD
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
=======
>>>>>>> 5b7361d93c9119cf1d2e9e6c885bed93f924d71b
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HOURLY ECONOMIC LOAD DISPATCH (Fitness Function)
%% Guided by : Dr. B.K.Panigrahi, V.Ravikumar Pandi, IIT Delhi
%% Coded by  : Krishnanand K.R., Santanu Kumar Nayak
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fitness Function to be called from the optimization algorithm
%% Evaluates a population of row vectors
function [y Count] = MY_FUNCTION(input_array)
siz = size(input_array,1);
Count = siz ;
y = zeros(siz,1);
for i =1:siz
    y(i,1) =  Fn_Eval(input_array(i,:));
end
end
%% Evaluates a single row vector
function y = Fn_Eval(x)
x=round(x*10000)/10000; %% For fixing the 4 digit precision
y = fn_ELD_6(x); %% Change Here to Call a Different System
end
%%
