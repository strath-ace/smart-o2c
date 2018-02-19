% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
clear all;
clc
format long e;
x =ones(30,2);
global initial_flag;
for i=1:28
   i 
   initial_flag = 0;
   [f,g,h]=CEC2017(x',i)
end
