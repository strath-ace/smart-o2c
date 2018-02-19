% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/. */
%
%------ Copyright (C) 2018 University of Strathclyde and Authors ------
%--------------- e-mail: smart@strath.ac.uk ---------------------------
%------------------- Authors: SMART developers team -------------------
function  [Cf,f]=fcon_red(fworst,C)
%
%  [Cf,f]=fcon_red(fworst,C)
%  
% augmented fitness function for constraint problems
%
%  INPUT
%            fin         : input fitness
%            fworst  :  reference fitness
%            C          : constraint vector
%
%  OUTPUT
%            f          : augmented fitness
%            Cf        : scalar function for constraint violation
%
% (c) Massimiliano Vasile 2003
%
Cf=0;
for i=1:length(C)
    Cf=Cf+max([0 C(i)]);
end
f=fworst+Cf;
return